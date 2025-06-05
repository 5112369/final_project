import sys
import os
import matplotlib.pyplot as plt
import unittest
import csv
import time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from strobealign import load_fasta, generate_strobemers, build_index, strobealign,plot_alignment_distribution,save_alignments_to_csv

class TestStrobeAlign(unittest.TestCase):

    def setUp(self):
        self.ref_path = "data/Ecoli_K12_MG1655.fa"
        self.query_path = "data/SRR5168216_1.fa"
        self.reference = load_fasta(self.ref_path)
        self.query = load_fasta(self.query_path)

    def test_save_alignments_to_csv(self):
        alignments = strobealign(self.query, self.reference, k=4, w_min=1, w_max=2, n=2, min_hits=1)
        print("比对结果数量:", len(alignments))

        # 生成带时间戳的文件名
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        ref_name = os.path.splitext(os.path.basename(self.ref_path))[0]
        query_name = os.path.splitext(os.path.basename(self.query_path))[0]
        out_csv = f"{ref_name}_vs_{query_name}_alignments.csv"
        out_png = f"{ref_name}_vs_{query_name}_distribution.png"

        save_alignments_to_csv(alignments, out_csv)
        print(f"比对结果已保存到: {out_csv}")

        plot_alignment_distribution(alignments, filename=out_png)
        print(f"分布图已保存到: {out_png}")

        self.assertGreater(len(alignments), 0, "没有找到任何比对结果，请检查输入数据或参数设置")
    def test_load_fasta(self):
        self.assertIsNotNone(self.reference)
        self.assertIsNotNone(self.query)

    def test_generate_strobemers(self):
        strobemers = generate_strobemers(self.query, k=5,w_min=2, w_max=4, n=2)
        print("strobemers数量:", len(strobemers))
        self.assertGreater(len(strobemers), 0)

    def test_build_index(self):
        index = build_index(self.reference, k=10)
        print("索引条目数:", len(index))
        self.assertGreater(len(index), 0)
    
if __name__ == '__main__':
       unittest.main()
    
    