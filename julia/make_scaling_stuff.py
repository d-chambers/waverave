"""
Simple python script to parse and display timing results
"""
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def parse_scaling_output(path):
    """Parse output into a dataframe"""
    with open(path, "r") as f:
        lines = f.readlines()
    data = []
    for line in lines:
        if line.startswith("Job"):
            continue
        grid_str = line.split("grid size of: (")[-1].split(")")[0]
        out = dict(
            nodes = int(line.split("node")[0].split("Used")[1].strip()),
            timesteps = int(line.split("time steps")[0].split("each with")[-1].strip()),
            grid_x = int(grid_str.split(",")[1].strip()),
            grid_y = int(grid_str.split(",")[0].strip()),
            time = float(line.split("time=")[1].split("seconds")[0].strip()),
        )
        data.append(out)
    df = pd.DataFrame(data)
    df['cells'] = df['grid_x'] * df['grid_y']
    return df



def plot_dataframe(df_weak, df_strong):
    """Make a plot of the input dataframes with timiing results"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))


    ax1.plot(df_strong['nodes'], df_strong['time'], 'o-')
    ax1.set_xlabel("Node Count")
    ax1.set_ylabel("Time (s)")
    ax1.set_title("Strong Scaling")
    ax1.grid(True)
    ax1.set_ylim(0, 1.1 * df_strong['time'].max())

    ax2.plot(df_weak['nodes'], df_weak['time'], 'o-')
    ax2.set_xlabel("Node Count")
    ax2.set_ylabel("Time (s)")
    ax2.set_title("Weak Scaling")
    ax2.grid(True)
    ax2.set_ylim(0, 1.1 * df_weak['time'].max())
    return fig, (ax1, ax2)


if __name__ == "__main__":
    strong_path = "strong_scaling.out"
    weak_path = "weak_scaling.out"
    report_path = Path(__file__).absolute().parent.parent / "writeup" / "images" / "julia"

    df_weak = parse_scaling_output(weak_path)
    df_strong = parse_scaling_output(strong_path)

    fig, *_ = plot_dataframe(df_weak, df_strong)
    fig.savefig(report_path / "scaling.png", dpi=300, 
                bbox_inches="tight", transparent=True)
    
    print("weak table")
    print(df_weak.to_markdown(index=False))
    print("strong table")
    print(df_strong.to_markdown(index=False))
