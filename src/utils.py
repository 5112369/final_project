from strobealign import load_fasta
def validate_fasta_format(filename):
    """Validate the format of a FASTA file."""
    with open(filename, 'r') as file:
        first_line = file.readline().strip()
        if not first_line.startswith('>'):
            raise ValueError(f"文件 {filename} 的格式不正确，应该以 '>' 开头。")
        
        for line in file:
            if not line.startswith('>') and not all(c in 'ACGTNacgtn' for c in line.strip()):
                raise ValueError(f"文件 {filename} 包含无效的字符。")
    
def read_fasta_file(filename):
    """Read a FASTA file and return the sequence."""
    validate_fasta_format(filename)
    sequence = load_fasta(filename)
    return sequence