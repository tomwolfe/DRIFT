import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

def create_dashboard(all_histories):
    """Creates a multi-panel Plotly dashboard from simulation results."""
    fig = make_subplots(
        rows=3, cols=1, 
        subplot_titles=(
            "Stochastic Signaling Dynamics (Normalized)", 
            "Metabolic Flux Outcome: Growth Rate (h⁻¹)",
            "Systemic Uncertainty Envelope (Phenotypic Response)"
        ),
        vertical_spacing=0.1
    )
    
    time = all_histories[0]['time']
    n_sims = len(all_histories)
    opacity = max(0.1, 1.0 / np.sqrt(n_sims)) if n_sims > 1 else 1.0
    
    # Trace Colors
    colors = {'PI3K': '#1f77b4', 'AKT': '#2ca02c', 'mTOR': '#d62728', 'Growth': '#ff7f0e'}

    # 1. Signaling Trajectories
    for i, hist in enumerate(all_histories):
        show_leg = (i == 0)
        # PI3K
        fig.add_trace(go.Scatter(
            x=time, y=hist['signaling'][:, 0], mode='lines', 
            line=dict(color=colors['PI3K']), opacity=opacity, 
            name='PI3K', legendgroup='sig', showlegend=show_leg
        ), row=1, col=1)
        # AKT
        fig.add_trace(go.Scatter(
            x=time, y=hist['signaling'][:, 1], mode='lines', 
            line=dict(color=colors['AKT']), opacity=opacity, 
            name='AKT', legendgroup='sig', showlegend=show_leg
        ), row=1, col=1)
        # mTOR
        fig.add_trace(go.Scatter(
            x=time, y=hist['signaling'][:, 2], mode='lines', 
            line=dict(color=colors['mTOR']), opacity=opacity, 
            name='mTOR', legendgroup='sig', showlegend=show_leg
        ), row=1, col=1)
        
        # 2. Growth Trajectories
        fig.add_trace(go.Scatter(
            x=time, y=hist['growth'], mode='lines', 
            line=dict(color=colors['Growth']), opacity=opacity, 
            name='Growth Rate', legendgroup='growth', showlegend=show_leg
        ), row=2, col=1)

    # 3. Uncertainty Envelope (Bottom Panel)
    growths = np.array([h['growth'] for h in all_histories])
    mean_growth = np.mean(growths, axis=0)
    std_growth = np.std(growths, axis=0)
    
    # Shade for std dev
    fig.add_trace(go.Scatter(
        x=np.concatenate([time, time[::-1]]),
        y=np.concatenate([mean_growth + std_growth, (mean_growth - std_growth)[::-1]]),
        fill='toself',
        fillcolor='rgba(255, 127, 14, 0.2)',
        line=dict(color='rgba(255,255,255,0)'),
        hoverinfo="skip",
        name='±1 Std Dev',
        showlegend=True
    ), row=3, col=1)

    # Mean line
    fig.add_trace(go.Scatter(
        x=time, y=mean_growth, mode='lines', 
        line=dict(color='black', width=2), 
        name='Mean Phenotype'
    ), row=3, col=1)

    # Update Axes
    fig.update_xaxes(title_text="Time (arbitrary units)", row=3, col=1)
    fig.update_yaxes(title_text="Normalized Activity", row=1, col=1)
    fig.update_yaxes(title_text="Objective Value", row=2, col=1)
    fig.update_yaxes(title_text="Ensemble Growth Rate", row=3, col=1)

    fig.update_layout(
        height=1000, 
        title_text="<b>DRIFT: Multi-Scale Stochastic Research Workbench</b><br>Drug-Target Interaction Propagation (PI3K/AKT/mTOR Axis)",
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    
    return fig
