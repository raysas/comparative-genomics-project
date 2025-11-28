#!/usr/bin/env python

import pandas as pd
from dash import Dash, html
import dash_bio as dashbio

# Load your CSV
df = pd.read_csv("../../output/statistics/duplicated_genes_info.csv")

# Prepare data for Circos
circos_data = []
for _, row in df.iterrows():
    circos_data.append({
        'block_id': str(row['chromosome']),    # chromosome segment
        'start': row['start_pos'],
        'end': row['end_pos'],
        'name': row['gene_id'],               # label
        'color': '#1f77b4'                    # any color
    })

# Prepare layout (chromosome sizes are approximated here)
chrom_sizes = df.groupby('chromosome')['end_pos'].max().to_dict()
layout = []
for chrom, size in chrom_sizes.items():
    layout.append({
        'id': str(chrom),
        'label': f'Chr{chrom}',
        'color': '#dddddd',
        'len': size
    })

# Initialize Dash app
app = Dash(__name__)

app.layout = html.Div([
    dashbio.Circos(
        id='genome-circos',
        layout=layout,
        tracks=[{
            'type': 'CHORDS',  # CHORDS if you have links, else 'HIGHLIGHT'
            'data': circos_data
        }],
        size=600
    )
])

if __name__ == '__main__':
    app.run(debug=True)

    