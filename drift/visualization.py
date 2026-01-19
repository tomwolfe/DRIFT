import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

def create_dashboard(results):
    """
    Creates a multi-panel Plotly dashboard from simulation results.

    Args:
        results (dict): Dictionary containing 'histories' and 'basal_growth'

    Returns:
        plotly.graph_objects.Figure: Interactive dashboard figure
    """
    all_histories = results.get('histories', [])
    basal_growth = results.get('basal_growth', 0.1) # Fallback

    if not all_histories:
        raise ValueError("all_histories cannot be empty")

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Stochastic Signaling Dynamics",
            "Drug Sensitivity: Kd vs Outcome",
            "Growth Rate Trajectories (%)",
            "Phenotypic Uncertainty Envelope (%)"
        ),
        vertical_spacing=0.15,
        horizontal_spacing=0.1,
        specs=[[{"type": "scatter"}, {"type": "scatter"}],
               [{"type": "scatter"}, {"type": "scatter"}]]
    )

    time = all_histories[0]['time']
    n_sims = len(all_histories)

    # 1. Sensitivity Analysis (Top Right)
    # Correlation between perturbed Kd and final growth
    perturbed_kds = [h['drug_kd'] for h in all_histories]
    final_growths_pct = [(h['growth'][-1] / basal_growth) * 100 for h in all_histories]

    fig.add_trace(go.Scatter(
        x=perturbed_kds, y=final_growths_pct,
        mode='markers',
        marker=dict(
            color=final_growths_pct,
            colorscale='Viridis',
            showscale=False,
            size=8,
            opacity=0.7,
            line=dict(width=1, color='DarkSlateGrey')
        ),
        name='Sensitivity Data',
        hovertemplate='Kd: %{x:.3f}<br>Growth: %{y:.1f}%<extra></extra>'
    ), row=1, col=2)

    # Add trendline for sensitivity
    z = np.polyfit(perturbed_kds, final_growths_pct, 1)
    p = np.poly1d(z)
    kd_range = np.linspace(min(perturbed_kds), max(perturbed_kds), 10)
    fig.add_trace(go.Scatter(
        x=kd_range, y=p(kd_range),
        mode='lines',
        line=dict(color='red', dash='dash'),
        name='Response Trend'
    ), row=1, col=2)

    # Trace Colors
    colors = {'PI3K': '#1f77b4', 'AKT': '#2ca02c', 'mTOR': '#d62728', 'Growth': '#ff7f0e'}

    # 2. Sampled Trajectories (Top Left & Bottom Left)
    max_traces = 3
    sample_indices = np.linspace(0, n_sims - 1, min(n_sims, max_traces), dtype=int)

    for i in sample_indices:
        hist = all_histories[i]
        show_leg = bool(i == sample_indices[0])
        
        # Signaling (Top Left)
        fig.add_trace(go.Scatter(
            x=time, y=hist['signaling'][:, 0], mode='lines',
            line=dict(color=colors['PI3K'], width=1.5), opacity=0.5,
            name='PI3K', legendgroup='pi3k', showlegend=show_leg
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=time, y=hist['signaling'][:, 1], mode='lines',
            line=dict(color=colors['AKT'], width=1.5), opacity=0.5,
            name='AKT', legendgroup='akt', showlegend=show_leg
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=time, y=hist['signaling'][:, 2], mode='lines',
            line=dict(color=colors['mTOR'], width=1.5), opacity=0.5,
            name='mTOR', legendgroup='mtor', showlegend=show_leg
        ), row=1, col=1)

        # Growth Trajectories normalized (Bottom Left)
        growth_pct = (hist['growth'] / basal_growth) * 100
        fig.add_trace(go.Scatter(
            x=time, y=growth_pct, mode='lines',
            line=dict(color=colors['Growth'], width=1), opacity=0.4,
            name='Growth (Sample)', legendgroup='growth', showlegend=show_leg
        ), row=2, col=1)

    # 3. Uncertainty Envelope (Bottom Right)
    growths_pct = np.array([(h['growth'] / basal_growth) * 100 for h in all_histories])
    mean_growth = np.mean(growths_pct, axis=0)
    std_growth = np.std(growths_pct, axis=0)

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
    ), row=2, col=2)

    # Mean line
    fig.add_trace(go.Scatter(
        x=time, y=mean_growth, mode='lines',
        line=dict(color='black', width=3),
        name='Mean Phenotype'
    ), row=2, col=2)

    # Metadata Annotation for User Context
    inhibition = all_histories[0].get('inhibition', 0) * 100
    mean_vitality = np.mean(final_growths_pct)

    fig.add_annotation(
        text=(f"<b>Summary Statistics</b><br>"
              f"Target Inhibition: {inhibition:.1f}%<br>"
              f"Mean Vitality: {mean_vitality:.1f}% of basal<br>"
              f"MC Iterations: {n_sims}"),
        xref="paper", yref="paper",
        x=1.15, y=0.5,
        showarrow=False,
        align="left",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="black",
        borderwidth=1,
        borderpad=10
    )

    # Update Axes
    # Help Annotation for New Users
    fig.add_annotation(
        text=(f"<b>How to Read This Dashboard</b><br>"
              f"• <b>Top Left:</b> Temporal signaling response. Proteins should stabilize.<br>"
              f"• <b>Top Right:</b> Sensitivity of growth to drug affinity (Kd).<br>"
              f"• <b>Bottom Left:</b> Individual cell trajectories showing stochastic drift.<br>"
              f"• <b>Bottom Right:</b> Statistical envelope of the population phenotype."),
        xref="paper", yref="paper",
        x=1.15, y=0.1,
        showarrow=False,
        align="left",
        bgcolor="rgba(230,240,255,0.9)",
        bordercolor="blue",
        borderwidth=1,
        borderpad=10
    )

    # Update Axes
    fig.update_xaxes(title_text="Time (Integration Steps)", row=1, col=1)
    fig.update_xaxes(title_text="Binding Affinity (Kd)", row=1, col=2)
    fig.update_xaxes(title_text="Time (Integration Steps)", row=2, col=1)
    fig.update_xaxes(title_text="Time (Integration Steps)", row=2, col=2)

    fig.update_yaxes(title_text="Protein Activity [0,1]", row=1, col=1)
    fig.update_yaxes(title_text="Final Growth (% Basal)", row=1, col=2)
    fig.update_yaxes(title_text="Growth Rate (% Basal)", row=2, col=1)
    fig.update_yaxes(title_text="Ensemble Growth (% Basal)", row=2, col=2)

    fig.update_layout(
        height=850,
        title_text="<b>DRIFT: Multi-Scale Stochastic Research Workbench</b><br>Multi-Scale Phenotypic Response Dashboard",
        template="plotly_white",
        legend=dict(orientation="v", yanchor="middle", y=0.5, xanchor="left", x=1.02),
        margin=dict(r=250, t=100) # Increased right margin for annotations
    )

    return fig
