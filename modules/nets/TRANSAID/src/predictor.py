#!/usr/bin/env python
# predictor.py

import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
from pathlib import Path
import multiprocessing as mp
from typing import List, Dict, Union, Tuple, Sequence, Optional
from dataclasses import dataclass
from Bio import SeqIO
from tqdm import tqdm
from functools import partial
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import torch.nn.functional as F

from .models import TRANSAID_Embedding
from .utils import sequence_to_tensor
from .orf_score import BayesianScorer


@dataclass
class PredictionResult:
    """预测结果的数据类"""

    sequence_id: str
    start_position: int
    stop_position: int
    tis_score: float
    tts_score: float
    kozak_score: float
    cai_score: float
    gc_score: float
    integrated_score: float
    protein_sequence: str
    nucleotide_sequence: str
    # 添加过滤状态字段
    passed_filter: bool = True
    filter_reason: str = ""


class FastaDataset(Dataset):
    """用于FASTA序列的数据集"""

    def __init__(self, fasta_path, max_len):
        self.records = list(SeqIO.parse(fasta_path, "fasta"))
        self.max_len = max_len

    def __len__(self):
        return len(self.records)

    def __getitem__(self, idx):
        record = self.records[idx]
        sequence = str(record.seq).upper().replace("U", "T")
        sequence_id = record.id

        # 转换为tensor
        tensor = sequence_to_tensor(sequence, self.max_len)
        return (
            tensor,
            sequence_id,
            min(len(sequence), self.max_len),
            sequence[: self.max_len],
        )


class TranslationPredictor:
    def __init__(
        self,
        model_path: str,
        device: Union[str, int] = "cpu",
        batch_size: int = 32,
        sequence_length: int = 27109,
        tis_cutoff: float = 0.1,
        tts_cutoff: float = 0.1,
        orf_length_cutoff: int = 50,
    ):
        """初始化预测器

        Args:
            model_path: 模型权重文件路径
            device: 计算设备,'cpu'或GPU设备ID(int)或'cuda'
            batch_size: 批处理大小
            sequence_length: 序列固定长度,默认27109
            tis_cutoff: 最小概率值,
            tts_cutoff: 最小概率值,
            orf_length_cutoff: 最小CDS长度
        """
        self.batch_size = batch_size
        self.sequence_length = sequence_length
        self.tis_cutoff = tis_cutoff
        self.tts_cutoff = tts_cutoff

        # 设置设备
        if isinstance(device, int) and device >= 0:
            self.device = f"cuda:{device}"
        elif device == "cuda" and torch.cuda.is_available():
            self.device = "cuda"
        else:
            self.device = "cpu"
        print("device", self.device)
        # 检查CUDA可用性
        if "cuda" in self.device and not torch.cuda.is_available():
            logging.warning("CUDA requested but not available. Falling back to CPU.")
            self.device = "cpu"
        print("device", self.device)
        logging.info(f"Using device: {self.device}")

        # 加载模型
        try:
            self.model = TRANSAID_Embedding()
            state_dict = torch.load(model_path, map_location=self.device)
            self.model.load_state_dict(state_dict)
            self.model.to(self.device)
            self.model.eval()
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {model_path}: {str(e)}")

        # 初始化遗传密码表
        self.genetic_code = {
            "ATA": "I",
            "ATC": "I",
            "ATT": "I",
            "ATG": "M",
            "ACA": "T",
            "ACC": "T",
            "ACG": "T",
            "ACT": "T",
            "AAC": "N",
            "AAT": "N",
            "AAA": "K",
            "AAG": "K",
            "AGC": "S",
            "AGT": "S",
            "AGA": "R",
            "AGG": "R",
            "CTA": "L",
            "CTC": "L",
            "CTG": "L",
            "CTT": "L",
            "CCA": "P",
            "CCC": "P",
            "CCG": "P",
            "CCT": "P",
            "CAC": "H",
            "CAT": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGA": "R",
            "CGC": "R",
            "CGG": "R",
            "CGT": "R",
            "GTA": "V",
            "GTC": "V",
            "GTG": "V",
            "GTT": "V",
            "GCA": "A",
            "GCC": "A",
            "GCG": "A",
            "GCT": "A",
            "GAC": "D",
            "GAT": "D",
            "GAA": "E",
            "GAG": "E",
            "GGA": "G",
            "GGC": "G",
            "GGG": "G",
            "GGT": "G",
            "TCA": "S",
            "TCC": "S",
            "TCG": "S",
            "TCT": "S",
            "TTC": "F",
            "TTT": "F",
            "TTA": "L",
            "TTG": "L",
            "TAC": "Y",
            "TAT": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGC": "C",
            "TGT": "C",
            "TGA": "*",
            "TGG": "W",
        }

        # 初始化评分器
        self.scorer = BayesianScorer(min_cds_len=orf_length_cutoff)

    def predict_batch(
        self,
        fasta_path: Union[str, Path],
        num_workers: int = 0,
        batch_size: int = None,
        pad_margin: int = 1024,
        max_tokens: int = 0,
    ) -> Sequence[PredictionResult]:
        """批量预测FASTA文件中的序列，使用真正的批处理

        Sequences are batched in length order and every batch is padded to its
        longest sequence plus ``pad_margin`` (rounded up to 256, capped at
        ``sequence_length``) instead of always to ``sequence_length``. The model
        is fully convolutional with a receptive field of 2,003 nt (at most 1,009
        positions to the right), so any margin >= 1,009 leaves the probabilities
        of every real position unchanged; results come back in input order as
        before. ``pad_margin < 0`` restores fixed-length padding in input order.
        ``max_tokens`` (> 0) caps batch size x padded length, mainly for GPUs.
        ``num_workers`` is kept for API compatibility; encoding runs in-process.
        """

        if batch_size is None:
            batch_size = self.batch_size
        dataset = FastaDataset(fasta_path, self.sequence_length)

        sequences = []
        for record in dataset.records:
            sequence = str(record.seq).upper().replace("U", "T")
            sequences.append((record.id, sequence[: self.sequence_length]))

        self.model.eval()
        per_sequence = [None] * len(sequences)  # (probs, transcript data, results)

        with (
            torch.no_grad(),
            tqdm(total=len(sequences), desc="Processing sequences") as pbar,
        ):
            for indices, padded_length in self._length_batches(
                [len(sequence) for _, sequence in sequences],
                batch_size,
                pad_margin,
                max_tokens,
            ):
                batch_tensors = torch.stack(
                    [sequence_to_tensor(sequences[k][1], padded_length) for k in indices]
                ).to(self.device)

                # 批量预测
                outputs = self.model(batch_tensors)  # [batch_size, seq_len, 3]
                probs = F.softmax(outputs, dim=-1)  # 使用 softmax 获取概率

                for i, k in enumerate(indices):
                    sequence_id, sequence = sequences[k]
                    per_sequence[k] = self._scan_sequence(
                        sequence_id, sequence, batch_tensors[i], probs[i]
                    )
                    pbar.update(1)

        all_results = []
        all_probs = []  # 收集所有预测的原始概率
        transcript_data = []  # 记录转录本ID和序列长度
        for seq_probs, data, results in per_sequence:
            all_probs.append(seq_probs)
            transcript_data.append(data)
            all_results.extend(results)

        # 按整合得分排序
        all_results.sort(key=lambda x: x.integrated_score, reverse=True)

        # 构建与训练阶段一致的结果格式
        formatted_probs = self._format_results_like_training(transcript_data, all_probs)

        return all_results, formatted_probs

    def _length_batches(self, lengths, batch_size, pad_margin, max_tokens):
        """Index batches and their padded length (see ``predict_batch``)."""
        if pad_margin < 0:
            return [
                (list(range(i, min(i + batch_size, len(lengths)))), self.sequence_length)
                for i in range(0, len(lengths), batch_size)
            ]

        def padded(length):
            return min(self.sequence_length, -(-(length + pad_margin) // 256) * 256)

        batches, current = [], []
        for k in sorted(range(len(lengths)), key=lambda k: lengths[k]):
            grown = current + [k]
            too_many = len(grown) > batch_size
            too_long = max_tokens > 0 and padded(lengths[k]) * len(grown) > max_tokens
            if current and (too_many or too_long):
                batches.append((current, padded(lengths[current[-1]])))
                grown = [k]
            current = grown
        if current:
            batches.append((current, padded(lengths[current[-1]])))
        return batches

    def _scan_sequence(self, sequence_id, sequence, tensor, probs):
        """Probabilities, transcript record and ORF results of one sequence."""
        effective_length = len(sequence)

        # 获取有效序列区域的预测概率
        valid_mask = tensor[:effective_length] != 0
        seq_probs = probs[:effective_length][valid_mask].cpu().numpy()

        data = {
            "transcript_id": sequence_id,
            "length": effective_length,
            "sequence": sequence,
        }

        # 提取预测概率
        tis_probs = probs[:effective_length, 0].cpu().numpy()
        tts_probs = probs[:effective_length, 1].cpu().numpy()

        # 寻找潜在的TIS和TTS位点
        potential_tis = []
        potential_tts = []

        # 寻找ATG和终止密码子
        for j in range(len(sequence) - 2):
            codon = sequence[j : j + 3]
            if codon == "ATG" and tis_probs[j : j + 3].mean() > self.tis_cutoff:
                potential_tis.append(j)
            elif (
                codon in ["TAA", "TAG", "TGA"]
                and tts_probs[j : j + 3].mean() > self.tts_cutoff
            ):
                potential_tts.append(j)

        # 这部分处理每个序列的ORF分析
        results = self._process_orfs(
            sequence_id,
            sequence,
            potential_tis,
            potential_tts,
            tis_probs,
            tts_probs,
        )
        return seq_probs, data, results

    def _process_orfs(
        self, sequence_id, sequence, potential_tis, potential_tts, tis_probs, tts_probs
    ):
        """处理ORFs并计算得分"""
        results = []
        for start in potential_tis:
            for stop in potential_tts:
                if stop > start and (stop - start) % 3 == 0:
                    cds = sequence[start : stop + 3]

                    # 计算各项得分
                    kozak_score = self.scorer.calculate_kozak_score(sequence, start)
                    cai_score = self.scorer.calculate_cai(cds)
                    gc_score = self.scorer.calculate_gc_score(cds)
                    tis_score = tis_probs[start : start + 3].mean()
                    tts_score = tts_probs[stop : stop + 3].mean()

                    # 计算整合得分
                    integrated_score = self.scorer.calculate_integrated_score(
                        tis_score, tts_score, kozak_score, cai_score, gc_score
                    )

                    # 翻译蛋白质序列
                    protein_seq = self._translate_sequence(cds)

                    result = PredictionResult(
                        sequence_id=sequence_id,
                        start_position=start,
                        stop_position=stop,
                        tis_score=tis_score,
                        tts_score=tts_score,
                        kozak_score=kozak_score,
                        cai_score=cai_score,
                        gc_score=gc_score,
                        integrated_score=integrated_score,
                        protein_sequence=protein_seq,
                        nucleotide_sequence=cds,
                    )
                    results.append(result)
        return results

    def _translate_sequence(self, sequence):
        """翻译核苷酸序列为蛋白质序列"""
        protein_seq = ""
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]
            if codon in self.genetic_code:
                protein_seq += self.genetic_code[codon]
        return protein_seq

    def _format_results_like_training(self, transcript_data, probabilities):
        """创建与训练阶段格式一致的结果"""
        results = []

        for i, trans_data in enumerate(transcript_data):
            probs = probabilities[i]
            # 创建兼容训练阶段格式的字典
            result = {
                "transcript_id": trans_data["transcript_id"],
                "predictions_probs": probs,  # 原始预测概率
                "length": trans_data["length"],
                # 其他可能需要的字段...
            }
            results.append(result)

        return results

    def to(self, device: Union[str, int]) -> None:
        """将模型移动到指定设备

        Args:
            device: 目标设备,'cpu'或GPU设备ID(int)或'cuda'
        """
        if isinstance(device, int) and device >= 0:
            target_device = f"cuda:{device}"
        elif device == "cuda" and torch.cuda.is_available():
            target_device = "cuda"
        else:
            target_device = "cpu"

        if target_device != self.device:
            self.device = target_device
            self.model.to(self.device)
            logging.info(f"Model moved to device: {self.device}")

    def get_device(self) -> str:
        """获取当前设备

        Returns:
            str: 当前设备名称
        """
        return self.device

    def is_gpu_available(self) -> bool:
        """检查是否有GPU可用

        Returns:
            bool: 是否有GPU可用
        """
        return torch.cuda.is_available()

    def get_available_devices(self) -> List[str]:
        """获取所有可用设备列表

        Returns:
            List[str]: 可用设备列表
        """
        devices = ["cpu"]
        if torch.cuda.is_available():
            devices.extend([f"cuda:{i}" for i in range(torch.cuda.device_count())])
        return devices

    def save_predictions(
        self,
        results: Sequence[PredictionResult],
        output_path: Union[str, Path],
        format: str = "csv",
    ) -> None:
        """保存预测结果

        Args:
            results: 预测结果列表
            output_path: 输出文件路径
            format: 输出格式,'csv'或'json'
        """
        output_path = Path(output_path)

        if format.lower() == "csv":
            import csv

            with open(output_path, "w", newline="") as f:
                writer = csv.writer(f)
                # 写入表头
                writer.writerow(
                    [
                        "sequence_id",
                        "start_position",
                        "stop_position",
                        "tis_score",
                        "tts_score",
                        "kozak_score",
                        "cai_score",
                        "gc_score",
                        "integrated_score",
                        "protein_sequence",
                        "nucleotide_sequence",
                    ]
                )
                # 写入数据
                for result in results:
                    writer.writerow(
                        [
                            result.sequence_id,
                            result.start_position,
                            result.stop_position,
                            result.tis_score,
                            result.tts_score,
                            result.kozak_score,
                            result.cai_score,
                            result.gc_score,
                            result.integrated_score,
                            result.protein_sequence,
                            result.nucleotide_sequence,
                        ]
                    )
        elif format.lower() == "json":
            import json

            with open(output_path, "w") as f:
                json.dump([vars(r) for r in results], f, indent=2)
        else:
            raise ValueError(f"Unsupported output format: {format}")

    # def __repr__(self) -> str:
    #    """返回预测器的字符串表示
    #
    #    Returns:
    #        str: 预测器的描述字符串
    #    """
    #    return (
    #        f"TranslationPredictor(device={self.device}, "
    #        f"batch_size={self.batch_size}, "
    #        f"sequence_length={self.sequence_length}, "
    #        f"probs_cutoff={self.probs_cutoff})"
    #    )
