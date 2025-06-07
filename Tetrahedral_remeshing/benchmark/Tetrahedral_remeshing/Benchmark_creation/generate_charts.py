import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import csv


def read_csv_with_auto_delimiter(csv_path):
    import csv as pycsv
    with open(csv_path, 'r', newline='', encoding='utf-8') as f:
        sample = f.read(4096)
        f.seek(0)
        sniffer = pycsv.Sniffer()
        try:
            dialect = sniffer.sniff(sample)
            delimiter = dialect.delimiter
        except pycsv.Error:
            delimiter = ';'  # Fallback to semicolon
        df = pd.read_csv(f, delimiter=delimiter)
    return df

def run(config, aggregated_data, csv_path=None):
    chart_files = []
    # Determine output_dir: default is 'Charts' subfolder next to CSV
    if csv_path is not None:
        csv_dir = os.path.dirname(os.path.abspath(csv_path))
        default_output_dir = os.path.join(csv_dir, 'Charts')
    else:
        default_output_dir = 'charts'
    output_dir = config.get('output_dir', default_output_dir)
    os.makedirs(output_dir, exist_ok=True)

    for chart in config['charts']:
        # Required fields
        for field in ['kind', 'x', 'y']:
            if field not in chart:
                raise ValueError(f"Each chart config must specify a '{field}' field.")
        kind = chart['kind']
        x = chart['x']
        y = chart['y']
        hue = chart.get('hue')
        title = chart.get('title', f"{kind.capitalize()} chart of {y} vs {x}")
        save_as = chart.get('save_as', f"{y}_vs_{x}_{kind}.png")
        save_path = os.path.join(output_dir, save_as)

        # --- Apply filter if present ---
        df = aggregated_data
        filter_str = ''
        if 'filter' in chart:
            filter_items = []
            for key, value in chart['filter'].items():
                df = df[df[key] == value]
                short_key = key.split('.')[-1]
                filter_items.append(f"{short_key}={value}")
            if filter_items:
                filter_str = " [filter: " + ", ".join(filter_items) + "]"
        if filter_str:
            title += filter_str

        plt.figure()
        if kind == 'scatter':
            sns.scatterplot(data=df, x=x, y=y, hue=hue)
        elif kind == 'box':
            sns.boxplot(data=df, x=x, y=y, hue=hue)
        elif kind == 'line':
            sns.lineplot(data=df, x=x, y=y, hue=hue)
        else:
            raise ValueError(f"Unsupported chart kind: {kind}")

        plt.title(title)
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()
        chart_files.append(save_path)
    return chart_files

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python generate_charts.py <config.yaml> <results.csv>")
        sys.exit(1)
    config_path = sys.argv[1]
    csv_path = sys.argv[2]
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    df = read_csv_with_auto_delimiter(csv_path)
    run(config, df, csv_path=csv_path) 