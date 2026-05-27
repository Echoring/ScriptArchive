# input format: TSV, only read 'Ks' header col
# mamba activate ks_analysis # pandas numpy matplotlib seaborn scipy sklearn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from sklearn.utils import resample
from sklearn.utils import check_random_state
from multiprocessing import Pool
import os
import argparse

def process_bootstrap_sample(args):
    ks_values, kdefactor, x_grid, peak_height, original_peaks_x, random_state = args
    rng = check_random_state(random_state)
    sample = resample(ks_values, random_state=rng)
    sample_kde = gaussian_kde(sample, bw_method=kdefactor)
    sample_peaks, _ = find_peaks(sample_kde(x_grid), height=peak_height*0.8)
    
    closest_peaks = []
    for peak_x in original_peaks_x:
        if len(sample_peaks) > 0:
            closest_idx = np.argmin(np.abs(x_grid[sample_peaks] - peak_x))
            closest_peaks.append(x_grid[sample_peaks[closest_idx]])
        else:
            closest_peaks.append(np.nan)
    return closest_peaks

def analyze_ks(file_path, bw_method, ks_cutoff, n_bootstrap, peak_height, prominence, distance, width):
    df = pd.read_csv(file_path, sep='\t', header=0)
    df['ks_numeric'] = pd.to_numeric(df['Ks'], errors='coerce')
    ks_values = df['ks_numeric'][(df['ks_numeric'] > 0) & (df['ks_numeric'] < ks_cutoff)].dropna()

    kde = gaussian_kde(ks_values, bw_method=bw_method)
    x_grid = np.linspace(0, ks_cutoff, 2000)
    kde_values = kde(x_grid)
    scale = ks_cutoff / (2000 - 1)
    peaks_indices, properties = find_peaks(kde_values, height=peak_height, prominence=prominence, distance=max(round(distance/scale), 1), width=max(round(width/scale), 1))
    original_peaks_x = x_grid[peaks_indices]
    
    peak_details = []
    if len(peaks_indices) > 0:
        with Pool(processes=os.cpu_count() - 1 or 1) as pool:
            args_list = [(ks_values, kde.factor, x_grid, peak_height, original_peaks_x, np.random.randint(1e6)) 
                        for _ in range(n_bootstrap)]
            bootstrap_results = pool.map(process_bootstrap_sample, args_list)  # shape: (n_bootstrap, n_peaks)
        
        bootstrap_peaks_all = np.array(bootstrap_results).T
        
        for idx, (peak_x, bootstrap_peaks) in enumerate(zip(original_peaks_x, bootstrap_peaks_all)):
            valid_peaks = bootstrap_peaks[~np.isnan(bootstrap_peaks)]
            bootstrap_std = np.std(valid_peaks) if len(valid_peaks) > 0 else np.nan
            
            peak_details.append({
                'peak_location': peak_x,
                'std_dev': bootstrap_std,
                'peak_height': properties['peak_heights'][idx],
                'peak_prominence': properties['prominences'][idx],
                'peak_width': properties['widths'][idx] * ks_cutoff/2000
            })

    return ks_values, x_grid, kde_values, peaks_indices, peak_details

def plot_ks_distributions(data, output_svg, xlim, peak_mark_cutoff):
    """
    Plots the Ks distributions and saves the figure.
    """
    plt.figure(figsize=(12, 8))
    colors = plt.colormaps['tab20'].resampled(len(data))

    for i, (genome, plot_data) in enumerate(data.items()):
        x_grid = plot_data['x_grid']
        kde_values = plot_data['kde_values']
        peaks = plot_data['peaks']

        sns.kdeplot(plot_data['ks_values'], label=genome, color=colors(i), bw_method=plot_data['bw_method'])

        if len(peaks) > 0:
            for peak_idx in peaks:
                peak_x = x_grid[peak_idx]
                if peak_x <= peak_mark_cutoff:
                    plt.axvline(x=peak_x, color=colors(i), linestyle='--', linewidth=1)
                    plt.text(peak_x + 0.005, kde_values[peak_idx], f'{peak_x:.3f}', fontsize=8, color=colors(i))

    plt.title("Ks Distribution")
    plt.xlabel('Ks Value')
    plt.ylabel('Density')
    plt.legend()
    plt.xlim(0, xlim)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(output_svg)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze and plot Ks distribution for paralog gene pairs.")
    parser.add_argument('--files', nargs='+', required=True, help='List of input Ks files, in the desired order for plotting.')
    parser.add_argument('--bw', default="scott", help='Bandwidth for Kernel Density Estimation: number / silverman / scott. (scott)')
    parser.add_argument('--distance', type=float, default=0, help='Minimum distance between Ks peaks. (0)')
    parser.add_argument('--width', type=float, default=0, help='Minimum width for peak detection. (0)')
    parser.add_argument('--height', type=float, default=1, help='Minimum height for peak detection. (1)')
    parser.add_argument('--prominence', type=float, default=0, help='Minimum prominence for peak detection. (0)')
    parser.add_argument('--ks_cutoff', type=float, default=1.0, help='Maximum Ks value to include in the analysis. (1)')
    parser.add_argument('--xlim', type=float, default=1.0, help='Upper limit for the x-axis of the plot. (1)')
    parser.add_argument('--peak_mark_cutoff', type=float, default=1.0, help='Ks value cutoff above which peak markers will not be shown. (1)')
    parser.add_argument('--bootstrap', type=int, default=1, help='Bootstrap to calculate std_dev. (1)')
    parser.add_argument('--output_svg', default='ks_distribution.svg', help='Output SVG file name. (ks_distribution.svg)')
    parser.add_argument('--output_csv', default='ks_analysis_summary.csv', help='Output CSV file name. (ks_analysis_summary.csv)')

    args = parser.parse_args()

    plot_data = {}
    summary_data = []

    for file in args.files:
        if not os.path.exists(file):
            print(f"Warning: File not found - {file}")
            continue

        genome_name = os.path.basename(file).split('.')[0]
        try:
            args.bw = float(args.bw)
        except:
            pass
        ks_values, x_grid, kde_values, peaks_indices, peak_details = analyze_ks(file, args.bw, args.ks_cutoff, args.bootstrap, args.height, args.prominence, args.distance, args.width)

        plot_data[genome_name] = {
            'ks_values': ks_values,
            'x_grid': x_grid,
            'kde_values': kde_values,
            'peaks': peaks_indices,
            'bw_method': args.bw,
            'peak_height': args.height
        }

        for peak in peak_details:
            summary_data.append({
                'Genome': genome_name,
                'Peak_Location': peak['peak_location'],
                'Peak_Std_Dev': peak['std_dev'],
                'Peak_Height': peak['peak_height'],
                'Peak_Prominence': peak['peak_prominence'],
                'Peak_Width': peak['peak_width'],
            })

    plot_ks_distributions(plot_data, args.output_svg, args.xlim, args.peak_mark_cutoff)

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(args.output_csv, index=False, float_format='%.5f')

    print(f"Analysis complete. Results saved to {args.output_svg} and {args.output_csv}")

if __name__ == '__main__':
    main()
