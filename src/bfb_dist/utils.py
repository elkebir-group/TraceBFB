CHROMOSOME_LENGTHS = {
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948, "chr4": 193574945,
    "chr5": 182045439, "chr6": 172126628, "chr7": 160567428, "chr8": 146259331,
    "chr9": 150617247, "chr10": 134758134, "chr11": 135127769, "chr12": 133324548,
    "chr13": 113566686, "chr14": 101161492, "chr15": 99753195, "chr16": 96330374,
    "chr17": 84276897, "chr18": 80542538, "chr19": 61707364, "chr20": 66210255,
    "chr21": 45090682, "chr22": 51324926, "chrX": 154259566, "chrY": 62460029,
    "chrM": 16569,
}

def format_red(s): return f"\033[91m{s}\033[00m"
def format_green(s): return f"\033[92m{s}\033[00m"
def format_yellow(s): return f"\033[93m{s}\033[00m"
def format_light_purple(s): return f"\033[94m{s}\033[00m"
def format_purple(s): return f"\033[95m{s}\033[00m"
def format_cyan(s): return f"\033[96m{s}\033[00m"
def format_light_gray(s): return f"\033[97m{s}\033[00m"
def format_black(s): return f"\033[90m{s}\033[00m"