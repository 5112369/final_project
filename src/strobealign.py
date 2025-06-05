from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import csv

def save_alignments_to_csv(alignments, filename="alignments.csv"):
    """将比对结果保存为CSV文件"""
    with open(filename, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["start of ref", "start of query ", "long"])
        for aln in alignments:
            writer.writerow(aln)
            
def plot_alignment_distribution(alignments, filename="alignment_distribution.png"):
    """绘制比对起点分布直方图，并添加详细批注"""
    if not alignments:
        print("没有比对结果，无法绘图。")
        return
    starts = [a[0] for a in alignments]
    counts, bins, patches = plt.hist(starts, bins=50, color='skyblue', edgecolor='black')
    plt.title("count of alignments by start position")
    plt.xlabel("count of alignments by start position")
    plt.ylabel("start position")
    max_count = max(counts)
    mean_start = sum(starts) / len(starts)
    min_start = min(starts)
    max_start = max(starts)
    plt.annotate(f"max_count: {int(max_count)}", xy=(0.7, 0.92), xycoords='axes fraction', fontsize=10, color='red')
    plt.annotate(f"mean_start: {mean_start:.2f}", xy=(0.7, 0.87), xycoords='axes fraction', fontsize=10, color='green')
    plt.annotate(f"min_start: {min_start}", xy=(0.7, 0.82), xycoords='axes fraction', fontsize=10, color='purple')
    plt.annotate(f"max_start: {max_start}", xy=(0.7, 0.77), xycoords='axes fraction', fontsize=10, color='orange')
    plt.annotate(f"alignments: {len(alignments)}", xy=(0.7, 0.72), xycoords='axes fraction', fontsize=10, color='blue')
    plt.axvline(mean_start, color='green', linestyle='dashed', linewidth=1)
    plt.text(mean_start, max_count/2, 'mean', color='green', rotation=90, va='center')
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()
def load_fasta(filename):
    
    """加载FASTA文件，返回第一个序列"""
    for record in SeqIO.parse(filename, "fasta"):
        return str(record.seq)
    return ""

def generate_strobemers(sequence, k=5, w_min=2, w_max=4, n=2):
    """生成strobemers"""
    if len(sequence) < k * n:
        return []
    
    strobemers = []
    total_length = len(sequence)
    
    for i in range(total_length - (k * n) + 1):
        strobes = []
        pos = i
        strobe = sequence[pos:pos+k]
        strobes.append((pos, strobe))
        
        for _ in range(n-1):
            window_start = pos + w_min
            window_end = min(pos + w_max, total_length - k)
            
            if window_start >= window_end:
                break
                
            next_pos = window_start + (hash(strobe) % (window_end - window_start))
            next_strobe = sequence[next_pos:next_pos+k]
            
            strobes.append((next_pos, next_strobe))
            pos = next_pos
            strobe = next_strobe
        
        if len(strobes) == n:
            strobemer_seq = ''.join([s[1] for s in strobes])
            strobemer_pos = tuple([s[0] for s in strobes])
            strobemers.append((strobemer_pos, strobemer_seq))
    
    return strobemers

def build_index(reference, k=15, w_min=20, w_max=30, n=3):
    """构建参考基因组索引"""
    index = defaultdict(list)
    strobemers = generate_strobemers(reference, k, w_min, w_max, n)
    
    for pos, strobemer in strobemers:
        index[strobemer].append(pos)
    
    return index

def strobealign(query, reference, k=15, w_min=20, w_max=30, n=3, min_hits=3):
    """执行比对"""
    ref_index = build_index(reference, k, w_min, w_max, n)
    query_strobemers = generate_strobemers(query, k, w_min, w_max, n)
    
    matches = defaultdict(list)
    for q_pos, q_strobemer in query_strobemers:
        if q_strobemer in ref_index:
            for r_pos in ref_index[q_strobemer]:
                offset = r_pos[0] - q_pos[0]
                matches[offset].append((r_pos, q_pos))
    
    results = []
    for offset, hit_list in matches.items():
        if len(hit_list) >= min_hits:
            r_start = min([pos[0][0] for pos in hit_list])
            q_start = min([pos[1][0] for pos in hit_list])
            r_end = max([pos[0][-1] + k for pos in hit_list])
            q_end = max([pos[1][-1] + k for pos in hit_list])
            
            match_length = min(r_end - r_start, q_end - q_start)
            results.append((r_start, q_start, match_length))
    
    results.sort(key=lambda x: -x[2])
    return results
