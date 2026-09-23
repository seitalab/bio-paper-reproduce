#!/usr/bin/env python3
"""Download all Figure 3 inputs: GEO processed matrices, IFN gene set, and ELISA source data."""
import argparse
import csv
import json
import zipfile
import xml.etree.ElementTree as ET
import gzip
from pathlib import Path
import shutil
import tarfile
import tempfile
from urllib.request import urlopen

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161336/suppl/GSE161336_RAW.tar"
PREFIXES = [
    "GSM4905024_mock_", "GSM4905025_d1_HSV-1_", "GSM4905026_d1_HSV-2_",
    "GSM4905027_d3_HSV-1_", "GSM4905028_d3_HSV-2_",
    "GSM4905029_d5_HSV-1_", "GSM4905030_d5_HSV-2_",
]
EXPECTED = {p + s for p in PREFIXES for s in
            ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")}


def check_gzip(path):
    # Reading to EOF verifies gzip CRC and detects truncated compressed data.
    with gzip.open(path, "rb") as stream:
        while stream.read(1024 * 1024):
            pass


def download_paper_sources(destination):
    """Keep original source files and extract the ELISA measurements with cell references."""
    destination.mkdir(parents=True, exist_ok=True)
    sources = {
        'elife-65762-fig3-data1-v1.xlsx': 'https://cdn.elifesciences.org/articles/65762/elife-65762-fig3-data1-v1.xlsx',
        'HALLMARK_INTERFERON_ALPHA_RESPONSE.mouse.json': 'https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName=HALLMARK_INTERFERON_ALPHA_RESPONSE&fileType=json',
    }
    for name, url in sources.items():
        path = destination / name
        if not path.exists():
            partial = path.with_suffix(path.suffix + '.part')
            with urlopen(url, timeout=120) as response, partial.open('wb') as output:
                shutil.copyfileobj(response, output)
            if name.endswith('.json'):
                json.loads(partial.read_text())['HALLMARK_INTERFERON_ALPHA_RESPONSE']['geneSymbols']
            else:
                with zipfile.ZipFile(partial) as workbook:
                    workbook.read('xl/worksheets/sheet1.xml')
            partial.replace(path)
        print(f'Paper source: {path}', flush=True)

    gene_set = json.loads((destination / 'HALLMARK_INTERFERON_ALPHA_RESPONSE.mouse.json').read_text())
    symbols = gene_set['HALLMARK_INTERFERON_ALPHA_RESPONSE']['geneSymbols']
    if not isinstance(symbols, list) or not symbols or not all(isinstance(g, str) and g for g in symbols):
        raise ValueError('Expected a nonempty list of mouse gene symbols in the MSigDB JSON')
    print(f'Validated IFN-alpha gene set: {len(symbols)} symbols (matching happens in the notebook).')

    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(destination / 'elife-65762-fig3-data1-v1.xlsx') as workbook:
        sheet = ET.fromstring(workbook.read('xl/worksheets/sheet1.xml'))
    cells = {cell.get('r'): cell.find('m:v', ns).text
             for cell in sheet.findall('.//m:c', ns)
             if cell.find('m:v', ns) is not None and cell.get('t') != 's'}
    path = destination / 'figure3E_ifnb_elisa.csv'
    with path.open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=['dpi', 'condition', 'ifnb_pg_ml', 'source_cell'])
        writer.writeheader()
        for day, start in [(4, 7), (5, 23)]:
            for column, condition in [('A', 'uninfected'), ('B', 'HSV-1'), ('C', 'HSV-2')]:
                for row in range(start, start + 9):
                    key = f'{column}{row}'
                    if key in cells:
                        writer.writerow({'dpi': day, 'condition': condition,
                                         'ifnb_pg_ml': float(cells[key]), 'source_cell': key})
    print(f'ELISA values extracted: {path}; author-assigned undetectable values of 1 retained.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data" / "ifn-neutrophil",
                        help="Root containing GSE161336/ and paper-source/ (default: project/data/ifn-neutrophil)")
    parser.add_argument("--data-dir", type=Path,
                        help="Directory for the archive and raw/ matrices")
    parser.add_argument("--paper-source-dir", type=Path,
                        help="Directory for Figure 3D gene set and Figure 3E ELISA data")
    args = parser.parse_args()
    data_root = args.data_root.expanduser().resolve()
    destination = (args.data_dir or data_root / 'GSE161336').expanduser().resolve()
    paper_source = (args.paper_source_dir or data_root / 'paper-source').expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / "GSE161336_RAW.tar"
    raw = destination / "raw"

    # 1. Download atomically; keep an existing archive and validate it below.
    if archive.exists():
        print(f"[1/3] Using existing archive: {archive}", flush=True)
    else:
        print(f"[1/3] Downloading {URL}", flush=True)
        partial = archive.with_suffix(".tar.part")
        try:
            with urlopen(URL, timeout=120) as response, partial.open("wb") as output:
                expected_size = response.headers.get("Content-Length")
                received = 0
                last_report = 0
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
                    received += len(chunk)
                    if received - last_report >= 10 * 1024 * 1024:
                        print(f"  Downloaded {received / 1024**2:.1f} MiB", flush=True)
                        last_report = received
            if expected_size is not None and received != int(expected_size):
                raise RuntimeError("Incomplete download: Content-Length mismatch")
            with tarfile.open(partial, "r:") as bundle:
                if {m.name for m in bundle.getmembers()} != EXPECTED:
                    raise RuntimeError("Unexpected archive contents; not promoting download")
            partial.replace(archive)
        except Exception:
            print(f"Download failed. Partial file: {partial}; rerunning restarts the download.", flush=True)
            raise

    # 2. Validate/extract only the expected regular files, without extractall.
    print("[2/3] Checking archive and extracting matrices", flush=True)
    raw.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:") as bundle:
        members = bundle.getmembers()
        if len(members) != len(EXPECTED) or {m.name for m in members} != EXPECTED:
            raise RuntimeError(f"Archive does not contain the expected 21 files: {archive}")
        with tempfile.TemporaryDirectory(prefix="extract-", dir=destination) as temp:
            staging = Path(temp)
            for member in members:
                if not member.isfile():
                    raise RuntimeError(f"Unexpected non-regular archive member: {member.name}")
                target = raw / member.name
                if target.is_file() and target.stat().st_size == member.size:
                    try:
                        check_gzip(target)
                        continue
                    except (OSError, EOFError):
                        pass
                staged = staging / member.name
                with bundle.extractfile(member) as source, staged.open("wb") as output:
                    shutil.copyfileobj(source, output)
                if staged.stat().st_size != member.size:
                    raise RuntimeError(f"Truncated file: {member.name}")
                check_gzip(staged)
                staged.replace(target)

    # 3. Download the additional Figure 3D/E inputs and report outputs.
    print("[3/3] Verified 21 gzip files across 7 samples.", flush=True)
    print(f"Archive retained: {archive}")
    print(f"Matrices: {raw}")
    download_paper_sources(paper_source)
    print("Python analysis uses all seven samples; paper Figure 3 caption specifies day 5.")


if __name__ == "__main__":
    main()
