import argparse
import nbformat
import subprocess
import os
from nbconvert import MarkdownExporter
from nbconvert.preprocessors import ExecutePreprocessor

def convert_percent_script_to_notebook(py_file, ipynb_file):
    subprocess.run(["pwd"], check=True)
    subprocess.run(["jupytext", "--to", "notebook", py_file, "--output", ipynb_file], check=True)

def execute_notebook(ipynb_file):
    with open(ipynb_file, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)
    
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": os.path.dirname(ipynb_file)}})
    
    with open(ipynb_file, "w", encoding="utf-8") as f:
        nbformat.write(nb, f)
        
from nbconvert.exporters import MarkdownExporter
from nbconvert.preprocessors import ExtractOutputPreprocessor

def export_notebook_to_markdown(ipynb_file, md_file, figures_dir=".figures"):
    import pathlib

    with open(ipynb_file, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    # Ensure figures_dir exists
    os.makedirs(figures_dir, exist_ok=True)

    # Set up exporter and preprocessor
    exporter = MarkdownExporter()
    exporter.output_files_dir = figures_dir

    # Use ExtractOutputPreprocessor to actually extract image outputs
    extract_preprocessor = ExtractOutputPreprocessor()
    extract_preprocessor.output_filename_template = os.path.join(figures_dir, "{unique_key}_{cell_index}_{index}{extension}")
    exporter.register_preprocessor(extract_preprocessor, enabled=True)

    # Export notebook
    (body, resources) = exporter.from_notebook_node(nb)

    # Save Markdown file
    with open(md_file, "w", encoding="utf-8") as f:
        f.write(body)

    # Save extracted resources (figures)
    for filename, data in resources.get("outputs", {}).items():
        if str(filename).startswith(figures_dir):
            print(filename)
            output_path = os.path.join(os.path.dirname(md_file), filename)
            pathlib.Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            pathlib.Path(output_path).write_bytes(data)


def main():
    import os
    parser = argparse.ArgumentParser(description="Convert a Python percent script to an executed notebook and Markdown.")
    parser.add_argument("input_py", help="Path to the input .py file (percent script)")
    parser.add_argument("output_md", help="Path to the output .md file")
    parser.add_argument("execute", action="store_true", help="Execute the notebook after conversion", default=False)
    args = parser.parse_args()

    args.input_py = os.path.abspath(args.input_py)
    args.output_md = os.path.abspath(args.output_md)

    ipynb_file = args.input_py.replace(".py", ".ipynb")

    if False:
        print(f"🔁 Converting {args.input_py} to notebook...")
        convert_percent_script_to_notebook(args.input_py, ipynb_file)

        print(f"⚙️  Executing notebook {ipynb_file}...")
        execute_notebook(ipynb_file)

    print(f"📝 Exporting to Markdown {args.output_md}...")
    export_notebook_to_markdown(ipynb_file, args.output_md)

    print("✅ Done!")

if __name__ == "__main__":
    main()