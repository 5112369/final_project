# Strobe Alignment Project

## Overview
The Strobe Alignment Project implements a strobe alignment algorithm for comparing biological sequences. This project provides tools to load FASTA files, generate strobemers, build an index for reference sequences, and perform alignments between query and reference sequences.

## Features
- Load FASTA files containing biological sequences.
- Generate strobemers for efficient sequence alignment.
- Build an index for quick lookups of reference sequences.
- Perform alignment and retrieve matching segments.

## Installation
To set up the project, ensure you have Python installed on your machine. Then, clone the repository and install the required dependencies:

```bash
git clone <repository-url>
cd strobealign-project
pip install -r requirements.txt
```

## Usage
To use the strobe alignment algorithm, you can run the main script located in the `src` directory. Ensure that your FASTA files are placed in the `data` directory.

Example command to run the alignment:

```bash
python src/strobealign.py
```

## Testing
Unit tests for the algorithm are provided in the `tests` directory. You can run the tests using a testing framework like `pytest`:

```bash
pytest tests/test_strobealign.py
```

## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue for any enhancements or bug fixes.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.