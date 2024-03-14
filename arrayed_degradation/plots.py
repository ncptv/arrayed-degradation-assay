import typing as tp
from io import BytesIO
import logging

import numpy as np
import numpy.typing as npt
import pandas as pd
import plotly.express as px  # type: ignore[import]
import plotly.graph_objects as go  # type: ignore[import]
import plotly.io as pio  # type: ignore[import]
from PIL import Image

pio.templates.default = "plotly_dark"
logging.getLogger("PIL").setLevel(logging.WARNING)  # removing PIL debug messages

PLOTLY_SCALE = 2
COLOR_OPTIONS = [
    "rgba(127, 200, 248, {opacity})",
    "rgba(204, 41, 54, {opacity})",
    "rgba(100, 201, 135, {opacity})",
    "rgba(255, 228, 94, {opacity})",
]


def plot_summary_chart(results: pd.DataFrame) -> Image.Image:
    fig = px.scatter(
        results,
        x="rna_id",
        y="half_life",
        error_y="half_life_std",
        title="Mean half-life",
        width=1200,
    )
    fig.update_layout(
        xaxis=dict(tickangle=-45, tickmode="linear"),
        yaxis=dict(title="Half-life in hours"),
        showlegend=False,
    )
    fig.update_traces(
        marker=dict(color=COLOR_OPTIONS[0].format(opacity=1), symbol="circle")
    )
    y_min = max(0, (results.half_life - results.half_life_std).min()) * 0.9
    y_max = min(20, (results.half_life + results.half_life_std).max()) * 1.1
    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor="gray", zeroline=True)
    fig.update_yaxes(
        showgrid=True,
        gridwidth=0.5,
        gridcolor="gray",
        zeroline=True,
        range=[y_min, y_max],
    )
    fig_bytes = fig.to_image(format="png", scale=PLOTLY_SCALE)
    return Image.open(BytesIO(fig_bytes))


class Plotter:
    """Plotter class for RNA stability assay."""

    def __init__(
        self,
        sequence_name: str,
        timepoints: list[float],
        n_replicate: int,
    ) -> None:
        self.sequence_name = sequence_name
        self.timepoint_to_color = {
            d: COLOR_OPTIONS[i % len(COLOR_OPTIONS)] for i, d in enumerate(timepoints)
        }
        self.max_epg_value = 0

        self.fig_decay_curve = go.Figure(
            layout=go.Layout(
                xaxis_title="Time (hours)",
                yaxis_title="Fraction remaining",
                showlegend=False,
                width=500,
            )
        )

        self.fig_residuals = go.Figure(
            layout=go.Layout(
                title="Exponent fit residuals",
                xaxis_title="Time (hours)",
                yaxis_title="Residual",
                width=350,
                yaxis_range=[-0.1, 0.1],
                showlegend=False,
            )
        )

        self.fig_areas = go.Figure(
            layout=go.Layout(
                title="RNA amount",
                xaxis_title="Time",
                yaxis_title="raw RNA amount (area under the peak)",
                width=350,
            )
        )

        # figure with electropherogram
        self.fig_epgs = [
            go.Figure(
                layout=go.Layout(
                    title=f"Sample EPG; replicate {i + 1}",
                    xaxis_title="Nucleotide length",
                    yaxis_title="Signal intensity",
                    xaxis_range=[0, 2500],
                    width=1200,
                )
            )
            for i in range(n_replicate)
        ]

    def plot_decay_curve(
        self,
        timepoints: npt.NDArray[np.float64],
        fractions_remaining: npt.NDArray[np.float64],
        fit_results: dict[str, float | tp.Callable],
        fit_fn: tp.Callable,
    ) -> None:
        self.fig_decay_curve.add_trace(
            go.Scatter(
                x=timepoints,
                y=fractions_remaining,
                name="data",
                mode="markers",
                marker=dict(color=COLOR_OPTIONS[2].format(opacity=1)),
            )
        )
        self.fig_decay_curve.add_trace(
            go.Scatter(
                x=np.linspace(min(timepoints), max(timepoints)),
                y=fit_fn(np.linspace(min(timepoints), max(timepoints))),
                name="fit",
                mode="lines",
                line=dict(color=COLOR_OPTIONS[0].format(opacity=1)),
            )
        )
        self.fig_decay_curve.update_layout(
            title=f"Decay curve of {self.sequence_name}<br>"
            f"half life {fit_results['half_life']:0.3f} h (+- {fit_results['half_life_std']:0.3f} std)<br>"
            f"R^2 {fit_results['r2_score']:0.2f}",
            yaxis_range=[0, 1.05],
        )

    def plot_residuals(
        self,
        timepoints: npt.NDArray[np.float64],
        fractions_remaining: npt.NDArray[np.float64],
        fit: tp.Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
    ) -> None:
        self.fig_residuals.add_trace(
            go.Scatter(
                x=timepoints,
                y=fractions_remaining - fit(timepoints),
                name="data",
                mode="markers",
                marker=dict(color=COLOR_OPTIONS[2].format(opacity=1)),
            )
        )
        self.fig_residuals.add_trace(
            go.Scatter(
                x=timepoints,
                y=np.zeros(len(timepoints)),
                name="fit",
                mode="lines",
                line=dict(color=COLOR_OPTIONS[0].format(opacity=1)),
            )
        )

    def plot_rna_amounts(
        self,
        timepoints: npt.NDArray[np.float64],
        control_mrna: npt.NDArray[np.float64] | None,
        peak_area: npt.NDArray[np.float64],
    ) -> None:
        if control_mrna is not None:
            self.fig_areas.add_trace(
                go.Scatter(
                    x=timepoints,
                    y=control_mrna,
                    name="control (m)RNA",
                    mode="markers",
                    line=dict(color="forestgreen"),
                )
            )
            self._add_line_with_std(
                self.fig_areas,
                timepoints,
                control_mrna,
                color="forestgreen",
            )
        self.fig_areas.add_trace(
            go.Scatter(
                x=timepoints,
                y=peak_area,
                name="peak area",
                mode="markers",
                line=dict(color=COLOR_OPTIONS[0].format(opacity=1)),
            )
        )
        self._add_line_with_std(
            self.fig_areas,
            timepoints,
            peak_area,
            color=COLOR_OPTIONS[0].format(opacity=1),
        )
        self.fig_areas.update_layout(
            legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99)
        )

    def plot_sample_trace(
        self,
        nucleotides: npt.NDArray[np.float64],
        trace: npt.NDArray[np.float64],
        timepoint: float,
        min_peak: float,
        max_peak: float,
        control_min_peak: float | None,
        control_max_peak: float | None,
        repl_no: int,
        remove_background: bool,
    ) -> None:
        # add epg
        self.fig_epgs[repl_no].add_trace(
            go.Scatter(
                x=nucleotides,
                y=trace,
                name=f"t={timepoint}h",
                mode="lines",
                line=dict(
                    width=3,
                    color=self.timepoint_to_color[timepoint].format(opacity=0.95),
                ),
            ),
        )
        # add target bounds
        self.fig_epgs[repl_no].add_vline(
            min_peak,
            line=dict(dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"),
        )
        self.fig_epgs[repl_no].add_vline(
            max_peak,
            line=dict(dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"),
        )
        # add area under the curve of target
        self.fig_epgs[repl_no].add_trace(
            self._get_area_trace(
                nucleotides=nucleotides,
                trace=trace,
                boundaries=(min_peak, max_peak),
                color=self.timepoint_to_color[timepoint].format(opacity=0.95),
                remove_background=remove_background,
            )
        )

        if control_min_peak is not None and timepoint == 0:
            # add control bounds
            self.fig_epgs[repl_no].add_vline(
                control_min_peak,
                line=dict(dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"),
            )
            self.fig_epgs[repl_no].add_vline(
                control_max_peak,
                line=dict(dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"),
            )
        if max(trace) > self.max_epg_value:
            self.max_epg_value = max(trace)
            for fig in self.fig_epgs:
                fig.update_yaxes(range=[0, self.max_epg_value])

    def save_fig(self) -> Image.Image:
        images = []
        for fig in [
            self.fig_decay_curve,
            self.fig_residuals,
            self.fig_areas,
            *self.fig_epgs,
        ]:
            buf = BytesIO()
            fig.write_image(buf, format="png", scale=PLOTLY_SCALE)
            buf.seek(0)
            images.append(Image.open(buf))
        decay_curve, residuals, areas = images.pop(0), images.pop(0), images.pop(0)
        big_image = self.concat_images([decay_curve, residuals, areas], axis=1)
        big_image = self.concat_images([big_image, *images], axis=0)
        return big_image

    def combined_plotting(
        self,
        fit_results: dict[str, float | tp.Callable],
        all_sample_timepoints: npt.NDArray[np.float64],
        all_control_areas: npt.NDArray[np.float64],
        all_peak_areas: npt.NDArray[np.float64],
        all_decay_normed_areas: npt.NDArray[np.float64],
    ) -> None:
        self.plot_rna_amounts(
            timepoints=all_sample_timepoints,
            control_mrna=all_control_areas if len(all_control_areas) != 0 else None,
            peak_area=all_peak_areas,
        )
        self.plot_decay_curve(
            timepoints=all_sample_timepoints,
            fractions_remaining=all_decay_normed_areas,
            fit_results=fit_results,
            fit_fn=tp.cast(tp.Callable, fit_results["fit"]),
        )
        self.plot_residuals(
            timepoints=all_sample_timepoints,
            fractions_remaining=all_decay_normed_areas,
            fit=tp.cast(tp.Callable, fit_results["fit"]),
        )

    @staticmethod
    def _get_area_trace(
        nucleotides: npt.NDArray[np.float64],
        trace: npt.NDArray[np.float64],
        boundaries: tuple[float, float],
        color: str,
        remove_background: bool,
    ) -> go.Scatter:
        mask = (nucleotides > boundaries[0]) & (nucleotides < boundaries[1])
        nucleotides = nucleotides[mask].copy()
        trace = trace[mask].copy()
        if remove_background:
            background_min = (np.abs(nucleotides - boundaries[0])).argmin()
            background_max = (np.abs(nucleotides - boundaries[1])).argmin()
            slope = (trace[background_max] - trace[background_min]) / (
                nucleotides[background_max] - nucleotides[background_min]
            )
            inter = trace[background_max] - slope * nucleotides[background_max]
            background_line = list(slope * nucleotides + inter)
        else:
            background_line = [0] * len(trace)
        y = background_line + trace[::-1].tolist()
        x = nucleotides.tolist() + nucleotides[::-1].tolist()
        return go.Scatter(
            x=x,
            y=y,
            fill="toself",
            mode="lines",
            line_color=color,
            line_width=0,
            showlegend=False,
        )

    @staticmethod
    def _add_line_with_std(
        fig: go.Figure, x: npt.NDArray, y: npt.NDArray, color: str
    ) -> None:
        tmp = pd.DataFrame({"x": x, "y": y})
        tmp = tmp.groupby("x")["y"].agg(["mean", "std"]).reset_index()
        fig.add_trace(
            go.Scatter(
                x=tmp["x"],
                y=tmp["mean"],
                mode="lines",
                line=dict(color=color),
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=tmp.x.tolist() + tmp.x.tolist()[::-1],
                y=(tmp["mean"] + tmp["std"]).tolist()
                + (tmp["mean"] - tmp["std"]).tolist()[::-1],
                fill="toself",
                mode="none",
                line_color=color,
                line_width=0,
                showlegend=False,
            )
        )

    @staticmethod
    def concat_images(images: list[Image.Image], axis: int) -> Image.Image:
        if axis == 0:  # Concatenate vertically
            widths, heights = zip(*(i.size for i in images))
            total_width = max(widths)
            total_height = sum(heights)
            new_image = Image.new("RGB", (total_width, total_height))
            y_offset = 0
            for image in images:
                new_image.paste(image, (0, y_offset))
                y_offset += image.size[1]
        elif axis == 1:  # Concatenate horizontally
            widths, heights = zip(*(i.size for i in images))
            total_width = sum(widths)
            total_height = max(heights)
            new_image = Image.new("RGB", (total_width, total_height))
            x_offset = 0
            for image in images:
                new_image.paste(image, (x_offset, 0))
                x_offset += image.size[0]
        return new_image
