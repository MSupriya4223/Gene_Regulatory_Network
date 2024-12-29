import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate
import plotly.graph_objs as go
from sklearn.preprocessing import MinMaxScaler
from dash import Dash, dcc, html, Input, Output, State, callback_context

# Mention the Dataset name
x = 'GSE18842.xlsx'

# Load the Dataset using pandas
df = pd.read_excel(x)
# print(df.head(5))

# Define the keys to specify 'Healthy' and 'Disease' samples
Healthy = 'Control'
Disease = 'Tumor'

# Preprocessing
df = df.drop(['!Sample_title'], axis=1)
df = df.dropna(subset=['Gene Symbol'])

# Handle the Duplicate rows based on Gene name
duplicates = df[df.duplicated(subset='Gene Symbol', keep=False)]
avg_duplicates = duplicates.groupby('Gene Symbol').mean().reset_index()
filtered_df = df[~df.duplicated(subset='Gene Symbol', keep=False)]
df = pd.concat([filtered_df, avg_duplicates])

# Normalize the matrix
scaler = MinMaxScaler()
normalized = scaler.fit_transform(df.iloc[:, 1:])
normalized_df = pd.DataFrame(normalized, columns=df.columns[1:]).reset_index(drop=True)
normalized_df['Gene Symbol'] = df['Gene Symbol'].reset_index(drop=True)
normalized_df.set_index('Gene Symbol', inplace=True)

# Transpose the matrix
transposed_df = normalized_df.transpose()

# 21 Genes that are obtained from venn diagram
top_selected_genes = ['SFTPC', 'DSG3', 'KRT14', 'KRT6B', 'MMP12', 'MMP1', 'GJB6', 'KRT5', 'SPRR1B', 'SPRR3', 'S100A2', 'AKR1B10', 'TMEM100', 'MT1M', 'SPRR1A', 'KRT16', 'WIF1', 'CLDN18', 'GJB2', 'JUP /// KRT17', 'KRT6A']

# Creating new data set with the top Genes
Disease_matrix = transposed_df[top_selected_genes]
Healthy_matrix = transposed_df[top_selected_genes]

Disease_matrix = Disease_matrix[Disease_matrix.index.str.contains(Disease)] 
Healthy_matrix = Healthy_matrix[Healthy_matrix.index.str.contains(Healthy)]

# Calculate correlations 
disease_corr = Disease_matrix.corr()
healthy_corr = Healthy_matrix.corr()

# print(disease_corr)
# print(healthy_corr)

# Set font scale for larger fonts
sns.set(font_scale=2.5)  # Adjust this value as needed for larger/smaller fonts

# Convert the correlation matrices to DataFrames
healthy_corr_df = pd.DataFrame(healthy_corr)
disease_corr_df = pd.DataFrame(disease_corr)

# Create and show heatmap for healthy samples with larger font
sns.clustermap(healthy_corr_df, cmap='coolwarm', method='ward', figsize=(16, 10),
               cbar_kws={'label': 'Correlation'})  # Colorbar label with larger font
               # xticklabels=1, yticklabels=1)  # Keep row and column labels visible
# plt.savefig('GSE18842_Healthy_Hierarchical_Clustering_Heatmap.png', dpi=500)
plt.show()

# Create and show heatmap for disease samples with larger font
sns.clustermap(disease_corr_df, cmap='coolwarm', method='ward', figsize=(16, 10),
               cbar_kws={'label': 'Correlation'})  # Colorbar label with larger font
               # xticklabels=1, yticklabels=1)  # Keep row and column labels visible
# plt.savefig('GSE18842_Disease_Hierarchical_Clustering_Heatmap.png', dpi=500)
plt.show()

# Create and visualize Disease network with Kamada-Kawai layout
G_disease = nx.Graph()

# Add nodes
for gene in top_selected_genes:
    G_disease.add_node(gene)

# Add edges with correlation values as weights
for i in range(len(top_selected_genes)):
    for j in range(i + 1, len(top_selected_genes)):
        corr_value = disease_corr.iloc[i, j]
        if abs(corr_value) > 0.3:  # Only add edges with significant correlation
            if corr_value > 0:
                if abs(corr_value) <= 0.5:
                    color = 'lightcoral'
                elif abs(corr_value) <= 0.7:
                    color = 'red'
                else:
                    color = 'indianred'
            else:
                if abs(corr_value) <= 0.5:
                    color = 'lightseagreen'
                elif abs(corr_value) <= 0.7:
                    color = 'blue'
                else:
                    color = 'blue'
            
            if abs(corr_value) <= 0.5:
                width = 1
            elif abs(corr_value) <= 0.7:
                width = 2
            else:
                width = 3
            
            G_disease.add_edge(top_selected_genes[i], top_selected_genes[j], weight=abs(corr_value), color=color, width=width)

# Compute Kamada-Kawai layout positions for Disease network
pos_disease = nx.kamada_kawai_layout(G_disease)

# Calculate node sizes based on degree (hub genes)
node_sizes_disease = [200 * G_disease.degree[node] for node in G_disease.nodes()]

# Draw the Disease network
edges_disease = G_disease.edges(data=True)
colors_disease = [edge[2]['color'] for edge in edges_disease]
widths_disease = [edge[2]['width'] for edge in edges_disease]

plt.figure(figsize=(16, 10))
nx.draw(G_disease, pos_disease, with_labels=True, node_color='wheat', edgecolors='black', node_size=node_sizes_disease, font_size=16)
# nx.draw_networkx_edge_labels(G_disease, pos_disease, edge_labels={(u, v): f'{d["weight"]:.2f}' for u, v, d in edges_disease})
nx.draw_networkx_edges(G_disease, pos_disease, edge_color=colors_disease, width=widths_disease)
plt.title('Disease Gene Correlation Network ')

# # Save the figure as a PNG file
# plt.savefig('GSE18842__Disease_GRN_Network.png', dpi=500)  # Adjust dpi for higher resolution if needed

# Calculate centrality measures for 'Disease' dataset
for i, gene1 in enumerate(disease_corr.index):
    for j, gene2 in enumerate(disease_corr.columns):
        if i != j:
            corr_value = disease_corr.iloc[i, j]
            if abs(corr_value) >= 0.3:  # Use absolute Threshold value as needed
                G_disease.add_edge(gene1, gene2, weight=corr_value)
                
# Calculate centrality measures
degree_centrality = nx.degree_centrality(G_disease)
closeness_centrality = nx.closeness_centrality(G_disease)
betweenness_centrality = nx.betweenness_centrality(G_disease)

disease_centrality_table = [[gene,
                        degree_centrality.get(gene, 0),
                        closeness_centrality.get(gene, 0),
                        betweenness_centrality.get(gene, 0)] for gene in disease_corr.index]

# Define header for the centrality measure table
headers = ["Gene", "Degree Centrality", "Closeness Centrality", "Betweenness Centrality"] 

print('\n\nCalculation of Centrality Measure: (Disease):\n')
print(tabulate(disease_centrality_table, headers=headers, tablefmt='grid', numalign='center'))

plt.show()

# Create and visualize Healthy network with Kamada-Kawai layout
G_healthy = nx.Graph()

# Add nodes
for gene in top_selected_genes:
    G_healthy.add_node(gene)

# Add edges with correlation values as weights
for i in range(len(top_selected_genes)):
    for j in range(i + 1, len(top_selected_genes)):
        corr_value = healthy_corr.iloc[i, j]
        if abs(corr_value) > 0.3:  # Only add edges with significant correlation
            if corr_value > 0:
                if abs(corr_value) <= 0.5:
                    color = 'lightgreen'
                elif abs(corr_value) <= 0.7:
                    color = 'green'
                else:
                    color = 'darkgreen'
            else:
                if abs(corr_value) <= 0.5:
                    color = 'lightseagreen'
                elif abs(corr_value) <= 0.7:
                    color = 'blue'
                else:
                    color = 'blue'
            
            if abs(corr_value) <= 0.5:
                width = 1
            elif abs(corr_value) <= 0.7:
                width = 2
            else:
                width = 3
            
            G_healthy.add_edge(top_selected_genes[i], top_selected_genes[j], weight=abs(corr_value), color=color, width=width)

# Compute Kamada-Kawai layout positions for Healthy network
pos_healthy = nx.kamada_kawai_layout(G_healthy)

# Calculate node sizes based on degree (hub genes)
node_sizes_healthy = [200 * G_healthy.degree[node] for node in G_healthy.nodes()]

# Draw the Healthy network
edges_healthy = G_healthy.edges(data=True)
colors_healthy = [edge[2]['color'] for edge in edges_healthy]
widths_healthy = [edge[2]['width'] for edge in edges_healthy]

plt.figure(figsize=(16, 10))
nx.draw(G_healthy, pos_healthy, with_labels=True, node_color='wheat', edgecolors='black', node_size=node_sizes_healthy, font_size=16)
# nx.draw_networkx_edge_labels(G_healthy, pos_healthy, edge_labels={(u, v): f'{d["weight"]:.2f}' for u, v, d in edges_healthy})
nx.draw_networkx_edges(G_healthy, pos_healthy, edge_color=colors_healthy, width=widths_healthy)
# plt.title('Healthy Gene Correlation Network')

# # Save the figure as a PNG file
# plt.savefig('GSE18842__Healthy_GRN_Network.png', dpi=500)  # Adjust dpi for higher resolution if needed

# Calculate centrality measures for 'Healthy' dataset
for i, gene1 in enumerate(healthy_corr.index):
    for j, gene2 in enumerate(healthy_corr.columns):
        if i != j:
            corr_value = healthy_corr.iloc[i, j]
            if abs(corr_value) >= 0.3:  # Use absolute Threshold value as needed
                G_healthy.add_edge(gene1, gene2, weight=corr_value)

# Calculate centrality measures
degree_centrality = nx.degree_centrality(G_healthy)
closeness_centrality = nx.closeness_centrality(G_healthy)
betweenness_centrality = nx.betweenness_centrality(G_healthy)

# Create a DataFrame for the centrality measures
healthy_centrality_df = pd.DataFrame({
    'Gene': healthy_corr.index,
    'Degree Centrality': [degree_centrality.get(gene, 0) for gene in healthy_corr.index],
    'Closeness Centrality': [closeness_centrality.get(gene, 0) for gene in healthy_corr.index],
    'Betweenness Centrality': [betweenness_centrality.get(gene, 0) for gene in healthy_corr.index]
})

# # Save the centrality measures to an Excel file
# healthy_centrality_df.to_excel('Healthy_Centrality_Measures.xlsx', index=False)

# Print the centrality measures table using tabulate
headers = ["Gene", "Degree Centrality", "Closeness Centrality", "Betweenness Centrality"]
print('\n\nCalculation of Centrality Measure: (Healthy):\n')
print(tabulate(healthy_centrality_df, headers=headers, tablefmt='grid', numalign='center'))

plt.show()

# Web host the network in ______spring layout_________
# This code can show the edges only whose node is clicked 
# And also show only weak, modrate amd stronger edges seperately 

grn_threshold = 0.3

# Create network graphs with reversed distances
def create_network_graph(corr_matrix, title, edge_filter=None, selected_node=None):
    G = nx.Graph()
    for i, gene1 in enumerate(corr_matrix.index):
        for j, gene2 in enumerate(corr_matrix.columns):
            if i != j:
                corr_value = corr_matrix.iloc[i, j]
                if abs(corr_value) >= grn_threshold:
                    G.add_edge(gene1, gene2, weight=corr_value)
    
    # Generate layout for the graph with spring layout algorithm using reversed distances
    pos = nx.spring_layout(G, weight=None, iterations=500, seed=42)
    
    edge_trace = []
    for edge in G.edges(data=True):
        corr_value = edge[2]['weight']
        if corr_value >= 0.75:
            color = 'indianred'
            width = 2
        elif corr_value >= 0.5:
            color = 'red'
            width = 1
        elif corr_value >= grn_threshold:
            color = 'lightcoral'
            width = 1
        elif corr_value <= -0.75:
            color = 'blue'
            width = 2
        elif corr_value <= -0.5:
            color = 'blue'
            width = 1
        elif corr_value <= -grn_threshold:
            color = 'lightseagreen'
            width = 0.5
        else:
            continue

        if edge_filter:
            if corr_value >= 0.75 and 'strong' not in edge_filter:
                continue
            if 0.5 <= corr_value < 0.75 and 'moderate' not in edge_filter:
                continue
            if grn_threshold <= corr_value < 0.5 and 'weak' not in edge_filter:
                continue
            if -0.75 <= corr_value < -0.5 and 'strong' not in edge_filter:
                continue
            if -0.5 <= corr_value < -0.75 and 'moderate' not in edge_filter:
                continue
            if -grn_threshold <= corr_value < -0.5 and 'weak' not in edge_filter:
                continue

        if selected_node and (edge[0] != selected_node and edge[1] != selected_node):
            continue
        
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace.append(go.Scatter(x=[x0, x1, None], y=[y0, y1, None],
                                     line=dict(width=width, color=color),
                                     hoverinfo='none',
                                     mode='lines'))
    
    node_trace = go.Scatter(x=[pos[k][0] for k in pos.keys()],
                            y=[pos[k][1] for k in pos.keys()],
                            text=[k for k in pos.keys()],
                            mode='markers+text',
                            textposition='top center',
                            hoverinfo='text',
                            # marker=dict(size=20, color='wheat'))
                            marker=dict(size=20, color='wheat',
                                        line=dict(width=.5, color='black')))  # Add border here
    
    fig = go.Figure(data=edge_trace + [node_trace],
                    layout=go.Layout(title=title, showlegend=False,
                                     hovermode='closest',
                                     margin=dict(b=20, l=5, r=5, t=40),
                                     xaxis=dict(showgrid=False, zeroline=False),
                                     yaxis=dict(showgrid=False, zeroline=False)))
    return fig

# Create heatmaps
def create_heatmap(corr_matrix, title, color_scheme):
    heatmap = go.Heatmap(z=corr_matrix.values,
                         x=corr_matrix.columns,
                         y=corr_matrix.index,
                         colorscale=color_scheme,
                         zmin=-1, zmax=1)
    layout = go.Layout(title=title, xaxis=dict(ticks=''), yaxis=dict(ticks=''))
    fig = go.Figure(data=[heatmap], layout=layout)
    return fig

# Initialize Dash app
app = Dash(__name__)

app.layout = html.Div([
    html.H1('Gene Network and Correlation Heatmap'),
    html.Div([
        html.Button('Toggle Weak Edges', id='weak-edges-btn', n_clicks=0),
        html.Button('Toggle Moderate Edges', id='moderate-edges-btn', n_clicks=0),
        html.Button('Toggle Strong Edges', id='strong-edges-btn', n_clicks=0),
    ]),
    dcc.Graph(id='disease-network', figure=create_network_graph(disease_corr, 'Disease Gene Network')),
    dcc.Graph(id='healthy-network', figure=create_network_graph(healthy_corr, 'Healthy Gene Network')),
    dcc.Graph(id='disease-heatmap', figure=create_heatmap(disease_corr, 'Disease Correlation Heatmap', 'RdBu')),
    dcc.Graph(id='healthy-heatmap', figure=create_heatmap(healthy_corr, 'Healthy Correlation Heatmap', 'RdBu')),
    dcc.Store(id='edge-filter', data=[]),
    dcc.Store(id='selected-node', data=None)
])

@app.callback(
    Output('edge-filter', 'data'),
    Input('weak-edges-btn', 'n_clicks'),
    Input('moderate-edges-btn', 'n_clicks'),
    Input('strong-edges-btn', 'n_clicks'),
    State('edge-filter', 'data')
)
def update_edge_filter(weak_clicks, moderate_clicks, strong_clicks, current_filter):
    ctx = callback_context
    if not ctx.triggered:
        return current_filter
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if button_id == 'weak-edges-btn':
        filter_type = 'weak'
    elif button_id == 'moderate-edges-btn':
        filter_type = 'moderate'
    elif button_id == 'strong-edges-btn':
        filter_type = 'strong'

    if filter_type in current_filter:
        current_filter.remove(filter_type)
    else:
        current_filter.append(filter_type)

    return current_filter

@app.callback(
    Output('disease-network', 'figure'),
    Output('healthy-network', 'figure'),
    Input('edge-filter', 'data'),
    Input('selected-node', 'data')
)
def update_network_graphs(edge_filter, selected_node):
    disease_fig = create_network_graph(disease_corr, ' ', edge_filter, selected_node)
    healthy_fig = create_network_graph(healthy_corr, ' ', edge_filter, selected_node)
    return disease_fig, healthy_fig

@app.callback(
    Output('selected-node', 'data'),
    Input('disease-network', 'clickData'),
    Input('healthy-network', 'clickData')
)
def update_selected_node(disease_click, healthy_click):
    ctx = callback_context
    if not ctx.triggered:
        return None
    graph_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if graph_id == 'disease-network':
        return disease_click['points'][0]['text'] if disease_click else None
    elif graph_id == 'healthy-network':
        return healthy_click['points'][0]['text'] if healthy_click else None
    return None

if __name__ == '__main__':
    app.run_server(debug=True)
