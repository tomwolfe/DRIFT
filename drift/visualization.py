import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np


def plot_death_diagnostics(history):
    """
    Visualizes the cause of metabolic collapse (Cell Death).
    Identifies which signaling-driven constraints were most restrictive.
    """
    if not history.get("cell_death", False):
        print("No cell death detected in this history. Diagnostic plot skipped.")
        return None

    death_step = history.get("death_step", 0)
    diag_msg = history.get("death_cause", "Unknown")
    
    # Extract signaling state at death
    signaling_at_death = history["signaling"][death_step]
    
    # We want to show which protein-reaction mapping was the 'killer'
    # This requires reaching into the solver's last known state or re-calculating
    # For the dashboard, we'll create a simple bar chart of constraints
    
    fig = go.Figure()
    
    # Plot signaling levels at time of death
    fig.add_trace(go.Bar(
        x=["PI3K", "AKT", "mTOR"], # Default species
        y=signaling_at_death,
        marker_color='red',
        name='Protein Level'
    ))

    fig.update_layout(
        title=f"<b>Cell Death Diagnostic</b><br>Cause: {diag_msg}",
        yaxis_title="Normalized Level [0,1]",
        xaxis_title="Signaling Species",
        template="plotly_white",
        shapes=[
            dict(
                type="line",
                yref="y", y0=0.2, y1=0.2,
                xref="paper", x0=0, x1=1,
                line=dict(color="Gray", dash="dash"),
            )
        ],
        annotations=[
            dict(
                x=0.5, y=0.1,
                text="Critical Threshold",
                showarrow=False,
                font=dict(color="Gray")
            )
        ]
    )
    
    return fig


def create_dashboard(results):
    """
    Creates a multi-panel Plotly dashboard from simulation results.

    Args:
        results (dict): Dictionary containing 'histories' and 'basal_growth'

    Returns:
        plotly.graph_objects.Figure: Interactive dashboard figure
    """
    all_histories = results.get("histories", [])
    basal_growth = results.get("basal_growth", 1.0) or 1.0  # Guard against zero/None

    if not all_histories:
        raise ValueError("all_histories cannot be empty")

    is_headless = all_histories[0].get("headless", False)

    # Pareto: Decimation for performance scaling
    original_time = all_histories[0]["time"]
    n_steps = len(original_time)
    max_display_points = 1000
    
    if n_steps > max_display_points:
        step_size = n_steps // max_display_points
        indices = np.arange(0, n_steps, step_size)
        time = original_time[indices]
        decimate = True
    else:
        time = original_time
        decimate = False

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "Stochastic Signaling Dynamics",
            "Drug Sensitivity: Kd vs Outcome",
            "Growth Rate Trajectories (%)",
            "Phenotypic Uncertainty Envelope (%)",
        ),
        vertical_spacing=0.15,
        horizontal_spacing=0.1,
        specs=[
            [{"type": "scatter"}, {"type": "scatter"}],
            [{"type": "scatter"}, {"type": "scatter"}],
        ],
    )

    n_sims = len(all_histories)

    # 1. Sensitivity Analysis (Top Right)
    # Correlation between perturbed Kd and final growth
    perturbed_kds = [h["drug_kd"] for h in all_histories]
    final_growths_pct = [(h["growth"][-1] / basal_growth) * 100 for h in all_histories]

    fig.add_trace(
        go.Scatter(
            x=perturbed_kds,
            y=final_growths_pct,
            mode="markers",
            marker=dict(
                color=final_growths_pct,
                colorscale="Viridis",
                showscale=False,
                size=8,
                opacity=0.7,
                line=dict(width=1, color="DarkSlateGrey"),
            ),
            name="Sensitivity Data",
            hovertemplate="Kd: %{x:.3f}<br>Growth: %{y:.1f}%<extra></extra>",
        ),
        row=1,
        col=2,
    )

    # Add trendline for sensitivity
    z = np.polyfit(perturbed_kds, final_growths_pct, 1)
    p = np.poly1d(z)
    kd_range = np.linspace(min(perturbed_kds), max(perturbed_kds), 10)
    fig.add_trace(
        go.Scatter(
            x=kd_range,
            y=p(kd_range),
            mode="lines",
            line=dict(color="red", dash="dash"),
            name="Response Trend",
        ),
        row=1,
        col=2,
    )

    # Trace Colors
    colors = {
        "PI3K": "#1f77b4",
        "AKT": "#2ca02c",
        "mTOR": "#d62728",
        "Growth": "#ff7f0e",
    }

    # 2. Sampled Trajectories (Top Left & Bottom Left)
    # Pareto: Only plot a limited number of traces to avoid Plotly lag
    max_traces = 5
    if n_sims > 50:
        sample_indices = np.linspace(0, n_sims - 1, max_traces, dtype=int)
    else:
        sample_indices = np.arange(min(n_sims, max_traces))

    for i in sample_indices:
        hist = all_histories[i]
        show_leg = bool(i == sample_indices[0])

        # Signaling (Top Left)
        for species_idx, (species_name, color) in enumerate(zip(["PI3K", "AKT", "mTOR"], [colors["PI3K"], colors["AKT"], colors["mTOR"]])):
            y_data = hist["signaling"][:, species_idx]
            if decimate:
                y_data = y_data[indices]
                
            fig.add_trace(
                go.Scatter(
                    x=time,
                    y=y_data,
                    mode="lines",
                    line=dict(color=color, width=1.5),
                    opacity=0.5,
                    name=species_name,
                    legendgroup=species_name.lower(),
                    showlegend=show_leg,
                ),
                row=1,
                col=1,
            )

        # Growth Trajectories normalized (Bottom Left)
        growth_pct = (hist["growth"] / basal_growth) * 100
        if decimate:
            growth_pct = growth_pct[indices]
            
        fig.add_trace(
            go.Scatter(
                x=time,
                y=growth_pct,
                mode="lines",
                line=dict(color=colors["Growth"], width=1),
                opacity=0.4,
                name="Growth (Sample)",
                legendgroup="growth",
                showlegend=show_leg,
            ),
            row=2,
            col=1,
        )

    # 3. Uncertainty Envelope (Bottom Right)
    # Memory-optimized stats calculation: Avoid one giant 2D matrix if n_sims is large
    if n_sims < 200:
        growths_pct_mat = np.array([(h["growth"] / basal_growth) * 100 for h in all_histories])
        mean_growth = np.mean(growths_pct_mat, axis=0)
        std_growth = np.std(growths_pct_mat, axis=0)
    else:
        # Incremental calculation to save memory
        sum_growth = np.zeros(n_steps)
        sum_growth_sq = np.zeros(n_steps)
        for h in all_histories:
            g = (h["growth"] / basal_growth) * 100
            sum_growth += g
            sum_growth_sq += g**2
        
        mean_growth = sum_growth / n_sims
        # Var = E[X^2] - (E[X])^2
        var_growth = (sum_growth_sq / n_sims) - (mean_growth**2)
        std_growth = np.sqrt(np.maximum(var_growth, 0)) # Max to avoid tiny numerical negatives
    
    if decimate:
        mean_growth = mean_growth[indices]
        std_growth = std_growth[indices]

    # Shade for std dev
    fig.add_trace(
        go.Scatter(
            x=np.concatenate([time, time[::-1]]),
            y=np.concatenate(
                [mean_growth + std_growth, (mean_growth - std_growth)[::-1]]
            ),
            fill="toself",
            fillcolor="rgba(255, 127, 14, 0.2)",
            line=dict(color="rgba(255,255,255,0)"),
            hoverinfo="skip",
            name="±1 Std Dev",
            showlegend=True,
        ),
        row=2,
        col=2,
    )

    # Mean line
    fig.add_trace(
        go.Scatter(
            x=time,
            y=mean_growth,
            mode="lines",
            line=dict(color="black", width=3),
            name="Mean Phenotype",
        ),
        row=2,
        col=2,
    )

    # Pareto: Add Basal Reference Line (100%) for visual context
    fig.add_hline(
        y=100.0,
        line_dash="dot",
        line_color="gray",
        annotation_text="Basal (100%)",
        annotation_position="bottom right",
        row=2,
        col=2,
    )

    # Metadata Annotation for User Context
    inhibition = all_histories[0].get("inhibition", 0) * 100
    mean_vitality = np.mean(final_growths_pct)
    
    status_msg = "METABOLISM: ACTIVE (COBRA)"
    status_color = "green"
    if is_headless:
        status_msg = "METABOLISM: HEADLESS (QUALITATIVE PROXY)"
        status_color = "red"

    fig.add_annotation(
        text=(
            f"<b>Summary Statistics</b><br>"
            f"Target Inhibition: {inhibition:.1f}%<br>"
            f"Mean Vitality: {mean_vitality:.1f}% of basal<br>"
            f"MC Iterations: {n_sims}<br>"
            f"<span style='color:{status_color}'><b>{status_msg}</b></span><br>"
            f"{'<i>Note: Proxy results not for publication.</i>' if is_headless else ''}"
        ),
        xref="paper",
        yref="paper",
        x=1.15,
        y=0.5,
        showarrow=False,
        align="left",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="black",
        borderwidth=1,
        borderpad=10,
    )

    # Update Axes
    # Help Annotation for New Users
    fig.add_annotation(
        text=(
            "<b>How to Read This Dashboard</b><br>"
            "• <b>Top Left:</b> Temporal signaling response. Proteins should stabilize.<br>"
            "• <b>Top Right:</b> Sensitivity of growth to drug affinity (Kd).<br>"
            "• <b>Bottom Left:</b> Individual cell trajectories showing stochastic drift.<br>"
            "• <b>Bottom Right:</b> Statistical envelope of the population phenotype."
        ),
        xref="paper",
        yref="paper",
        x=1.15,
        y=0.1,
        showarrow=False,
        align="left",
        bgcolor="rgba(230,240,255,0.9)",
        bordercolor="blue",
        borderwidth=1,
        borderpad=10,
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

    title_suffix = " (HEADLESS - QUALITATIVE PROXY)" if is_headless else ""
    fig.update_layout(
        height=850,
        title_text=f"<b>DRIFT: Multi-Scale Stochastic Research Workbench{title_suffix}</b><br>Multi-Scale Phenotypic Response Dashboard",
        template="plotly_white",
        legend=dict(orientation="v", yanchor="middle", y=0.5, xanchor="left", x=1.02),
        margin=dict(r=250, t=100),  # Increased right margin for annotations
    )

    return fig

