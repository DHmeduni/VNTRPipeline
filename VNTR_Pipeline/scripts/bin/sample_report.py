#!/usr/bin/env python3
import sys
import subprocess
import importlib
import os
import re


# --- Dependency check and install ---
def ensure_package(pkg_name, import_name=None, min_version=None, max_version=None, reinstall_if_bad=False):
    import_name = import_name or pkg_name
    from packaging import version as pv
    try:
        pkg = importlib.import_module(import_name)
        v = getattr(pkg, "__version__", "0")
        if (max_version and pv.parse(v) >= pv.parse(max_version)) or (min_version and pv.parse(v) < pv.parse(min_version)):
            if reinstall_if_bad:
                print(f"⚠️  {pkg_name} {v} incompatible — reinstalling...")
                subprocess.check_call([sys.executable, "-m", "pip", "uninstall", "-y", pkg_name])
                spec = pkg_name
                if max_version:
                    spec += f"<{max_version}"
                if min_version:
                    spec += f">={min_version}"
                subprocess.check_call([sys.executable, "-m", "pip", "install", spec])
                pkg = importlib.import_module(import_name)
        return pkg
    except ImportError:
        spec = pkg_name
        if max_version:
            spec += f"<{max_version}"
        if min_version:
            spec += f">={min_version}"
        print(f"📦 Installing missing package: {spec}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", spec])
        return importlib.import_module(import_name)

# Ensure required packages
pandas = ensure_package("pandas")
openpyxl = ensure_package("openpyxl")
reportlab = ensure_package("reportlab", max_version="3.6", reinstall_if_bad=True)
biopython = ensure_package("biopython", import_name="Bio")


# Now safe to import the rest of ReportLab
from Bio.Seq import Seq
from datetime import datetime
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import KeepTogether

def dna_to_protein(seq):
    try:
        return str(Seq(seq).translate(to_stop=True))
    except:
        return "INVALID"

def truncate_sample_name(name):
    #parts = re.split(r"[_-]", name)
    #if len(parts) >= 3:
        #return "".join(parts[:3])   # part1 + sep1 + part2
    #else:
        #return name
    pattern = r"^([^_-]+)([_-])([^_-]+)"   # first part, separator, second part
    match = re.match(pattern, name)
    if match:
        return match.group(1) + match.group(2) + match.group(3)
    return name

def safe_var(var_name, default="N/A"):
    return locals().get(var_name, default) or default

def main():
    if len(sys.argv) != 3:
        print("Usage: python sample_report.py <sample_dir> <sample_sheet.csv>")
        sys.exit(1)

    sample_dir = os.path.abspath(sys.argv[1])
    raw_sample_name = os.path.basename(sample_dir)
    sample_name = truncate_sample_name(raw_sample_name)
    sample_sheet = os.path.abspath(sys.argv[2])

    if not os.path.exists(sample_sheet):
        print(f"Error: sample_sheet.csv not found at {sample_sheet}")

    # --- Read CSV (tab-delimited) ---
    if os.path.exists(sample_sheet):
        df_sheet = pandas.read_csv(sample_sheet, sep='\t', header=None)
        operator = experiment_name = experiment_date = ''
        if len(df_sheet) > 1:
            operator = str(df_sheet.iloc[1, 0]) if not pandas.isna(df_sheet.iloc[1, 0]) else ''
            experiment_name = str(df_sheet.iloc[1, 1]) if df_sheet.shape[1] > 1 and not pandas.isna(df_sheet.iloc[1, 1]) else ''
            experiment_date = str(df_sheet.iloc[1, 2]) if df_sheet.shape[1] > 2 and not pandas.isna(df_sheet.iloc[1, 2]) else ''

    # Lookup barcode using truncated sample name
    
        df_data = pandas.read_csv(sample_sheet, sep='\t')
        name_col = 'Sample_Name' if 'Sample_Name' in df_data.columns else df_data.columns[0]
        sample_row = df_data[df_data[name_col] == sample_name]
        barcode = str(sample_row.iloc[0, 5]) if not sample_row.empty and sample_row.shape[1] > 5 and not pandas.isna(sample_row.iloc[0, 5]) else sample_name

    # --- Determine Result based on new_and_lof_seqs.xlsx ---
    result_text = "UNREMARKABLE"
    lof_file = os.path.join(sample_dir, "new_and_lof_seqs.xlsx")
    if os.path.exists(lof_file):
        try:
            df_lof = pandas.read_excel(lof_file)
            if df_lof.shape[1] >= 2 and any(df_lof.iloc[:, 1].astype(str).str.contains("L", na=False)):
                result_text = "FS/Nonsense detected"
        except Exception as e:
            print(f"⚠️ Could not process {lof_file}: {e}")

    # --- Prepare PDF ---
    pdf_path = os.path.join(sample_dir, f"{sample_name}_report.pdf")
    doc = SimpleDocTemplate(pdf_path, pagesize=A4)
    styles = getSampleStyleSheet()
    elements = []


    today = datetime.today().strftime("%Y-%m-%d")

    # dictionary with header info
    data = {
    "sample_name": locals().get('sample_name', 'N/A'),
    "experiment_name": locals().get('experiment_name', 'N/A'),
    "experiment_date": locals().get('experiment_date', 'N/A'),
    "operator": locals().get('operator', 'N/A'),
    "barcode": locals().get('barcode', 'N/A'),
    "today": locals().get('today', 'N/A')
    }
    
    # Header
    elements.append(Paragraph("<b>MUC1 VNTR PCR Report</b>", styles["Title"]))
    elements.append(Spacer(1, 12))
    elements.append(Paragraph(f"<b>Sample Name:</b> {data['sample_name']}", styles["Normal"]))
    elements.append(Paragraph(f"<b>Experiment:</b> {data['experiment_name']}", styles["Normal"]))
    elements.append(Paragraph(f"<b>Experiment Date:</b> {data['experiment_date']}", styles["Normal"]))
    elements.append(Paragraph(f"<b>Operator:</b> {data['operator']}", styles["Normal"]))
    elements.append(Paragraph(f"<b>Barcode:</b> {data['barcode']}", styles["Normal"]))
    elements.append(Paragraph(f"<b>Report Generated:</b> {data['today']}", styles["Normal"]))
    # --- Color the result text if LOSS OF FUNCTION ---
    if result_text == "FS/Nonsense detected":
        colored_result = f"<b>Result:</b> <font color='red'>{result_text or 'N/A'}</font>"
    else:
        colored_result = f"<b>Result:</b> {result_text or 'N/A'}"

    elements.append(Paragraph(colored_result, styles["Normal"]))
    elements.append(Spacer(1, 24))

    for excel_name in ["vntr_length.xlsx"]:
        excel_path = os.path.join(sample_dir, excel_name)
        if os.path.exists(excel_path):
            df = pandas.read_excel(excel_path)
            table_data = [list(df.columns)] + df.values.tolist()
            table = Table(table_data, repeatRows=1)
            table.setStyle(TableStyle([
                ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
                ('GRID', (0,0), (-1,-1), 0.25, colors.black),
                ('FONTSIZE', (0,0), (-1,-1), 8),
            ]))
            section = [
                Paragraph(f"<b>{excel_name}</b>", styles["Heading3"]),
                table,
                Spacer(1, 24)
            ]
            elements.append(KeepTogether(section))

    # Add images
    for img_name in ["best_trviz_fig.png", "best_trviz_protein_fig.png", None]:
        if img_name is None:
            for f in os.listdir(sample_dir):
                if f.endswith("histogram.png"):
                    img_name = f
                    break
        if img_name:
            img_path = os.path.join(sample_dir, img_name)
            if os.path.exists(img_path):
                section = [
                    Paragraph(f"<b>{img_name}</b>", styles["Heading3"]),
                    Image(img_path, width=600, height=350),
                    Spacer(1, 24)
                ]
                elements.append(KeepTogether(section))

    # Add Excel tables
    for excel_name in ["new_and_lof_seqs.xlsx", "seqs_distribution.xlsx"]:
        excel_path = os.path.join(sample_dir, excel_name)

        if os.path.exists(excel_path):

            df = pandas.read_excel(excel_path)  # <-- or pandas.read_excel if you kept that name

            if excel_name == "new_and_lof_seqs.xlsx":

                # Column 4 = index 3
                protein_df = df.iloc[:, [0, 1, 2]].copy()
                protein_df["Protein sequence"] = df.iloc[:, 3].astype(str).apply(dna_to_protein)

                # ---------------- DNA table ----------------
                table_data = [list(df.columns)] + df.values.tolist()
                table = Table(table_data, repeatRows=1)

                table.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
                    ('GRID', (0,0), (-1,-1), 0.25, colors.black),
                    ('FONTSIZE', (0,0), (-1,-1), 8),
                ]))

                elements.append(KeepTogether([
                    Paragraph(f"<b>{excel_name}</b>", styles["Heading3"]),
                    table,
                    Spacer(1, 24)
                ]))


                # ---------------- Protein table ----------------

                protein_table_data = [list(protein_df.columns)] + protein_df.values.tolist()

                protein_table = Table(protein_table_data, repeatRows=1)

                protein_table.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
                    ('GRID', (0,0), (-1,-1), 0.25, colors.black),
                    ('FONTSIZE', (0,0), (-1,-1), 8),
                ]))

                elements.append(KeepTogether([
                    Paragraph("<b>Protein sequences</b>", styles["Heading3"]),
                    protein_table,
                    Spacer(1, 24)
                ]))


            else:
                # ---------------- Other Excel files normal ----------------
                table_data = [list(df.columns)] + df.values.tolist()
                table = Table(table_data, repeatRows=1)

                table.setStyle(TableStyle([
                    ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
                    ('GRID', (0,0), (-1,-1), 0.25, colors.black),
                    ('FONTSIZE', (0,0), (-1,-1), 8),
                ]))

                elements.append(KeepTogether([
                    Paragraph(f"<b>{excel_name}</b>", styles["Heading3"]),
                    table,
                    Spacer(1, 24)
                ]))


    doc.build(elements)
    print(f"✅ PDF created: {pdf_path}")

if __name__ == "__main__":
    main()




