import logging
import typing as tp
from io import BytesIO

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
        rna_data: pd.DataFrame,
        fit_results: dict[str, float | tp.Callable],
        remove_background: bool,
        control_peak_min: float | None = None,
        control_peak_max: float | None = None,
    ) -> None:
        self.rna_data = rna_data
        self.fit_results = fit_results
        self.remove_background = remove_background
        self.control_peak_min = control_peak_min
        self.control_peak_max = control_peak_max

        self.fit_fn: tp.Callable = tp.cast(tp.Callable, fit_results["fit"])
        self.sequence_name = rna_data.rna_id.unique()[0]
        self.timepoints = sorted(rna_data.timepoint.unique())
        self.timepoint_to_color = {
            d: COLOR_OPTIONS[i % len(COLOR_OPTIONS)]
            for i, d in enumerate(self.timepoints)
        }
        self.replicates = sorted(rna_data.replicate.unique())
        self.max_epg_value = 0

        self.fig_decay_curve = go.Figure(
            layout=go.Layout(
                xaxis_title="Time (hours)",
                yaxis_title="Fraction remaining",
                width=500,
                legend=dict(
                    yanchor="top",
                    y=1.2,
                    xanchor="right",
                    x=1.15,
                    font_family="Courier New",
                    font_size=10,
                ),
            )
        )

        self.fig_residuals = go.Figure(
            layout=go.Layout(
                title="Exponent fit residuals",
                xaxis_title="Time (hours)",
                yaxis_title="Residual",
                width=350,
                showlegend=False,
            )
        )

        self.fig_areas = go.Figure(
            layout=go.Layout(
                title="RNA amount",
                xaxis_title="Time",
                yaxis_title="RNA amount (area under the peak)",
                width=350,
            )
        )

        self.fig_epgs = [
            go.Figure(
                layout=go.Layout(
                    title=f"Sample EPG; replicate {i + 1}",
                    xaxis_title="Nucleotide length",
                    yaxis_title="Signal intensity",
                    xaxis_range=[0, 2500],
                    width=1200,
                    legend=dict(
                        yanchor="bottom",
                        y=1,
                        xanchor="left",
                        x=0.5,
                        font_family="Courier New",
                        font_size=10,
                    ),
                )
            )
            for i in range(len(self.replicates))
        ]

    def plot_decay_curve(self) -> None:
        for (repl, timepoint), data in self.rna_data.groupby(
            ["replicate", "timepoint"]
        ):
            self.fig_decay_curve.add_trace(
                go.Scatter(
                    x=data.timepoint,
                    y=data.fraction_remaining,
                    mode="markers",
                    name=f"repl {repl+1}; tp {timepoint: <5}h",
                    marker=dict(
                        color=self.timepoint_to_color[timepoint].format(opacity=1),
                        symbol=repl,
                    ),
                )
            )
        self.fig_decay_curve.add_trace(
            go.Scatter(
                x=np.linspace(min(self.timepoints), max(self.timepoints)),
                y=self.fit_fn(np.linspace(min(self.timepoints), max(self.timepoints))),
                name="fit",
                mode="lines",
                line=dict(color=COLOR_OPTIONS[0].format(opacity=1)),
            )
        )
        self.fig_decay_curve.update_layout(
            title=f"Decay curve of {self.sequence_name}<br>"
            f"half life {self.fit_results['half_life']:0.3f} h (+- {self.fit_results['half_life_std']:0.3f} std)<br>"
            f"R^2 {self.fit_results['r2_score']:0.2f}",
            yaxis_range=[0, 1.05],
        )

    def plot_residuals(self) -> None:
        timepoints = self.rna_data.timepoint

        self.fig_residuals.add_trace(
            go.Scatter(
                x=timepoints,
                y=self.rna_data.fraction_remaining - self.fit_fn(timepoints),
                name="data",
                mode="markers",
                marker=dict(
                    color=[
                        self.timepoint_to_color[i].format(opacity=1)
                        for i in self.rna_data.timepoint
                    ],
                    symbol=self.rna_data.replicate,
                ),
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

    def plot_rna_amounts(self) -> None:
        timepoints = self.rna_data.timepoint
        peak_area = self.rna_data.target_peak_areas
        self.fig_areas.add_trace(
            go.Scatter(
                x=timepoints,
                y=peak_area,
                name="peak area",
                mode="markers",
                marker=dict(
                    color=[
                        self.timepoint_to_color[i].format(opacity=1)
                        for i in self.rna_data.timepoint
                    ],
                    symbol=self.rna_data.replicate,
                ),
            )
        )
        self._add_line_with_std(
            self.fig_areas,
            timepoints,
            peak_area,
            color=COLOR_OPTIONS[0].format(opacity=1),
        )

    def plot_sample_traces(self) -> None:
        for (repl, timepoint), data in self.rna_data.groupby(
            ["replicate", "timepoint"]
        ):
            if "epg_normed" in data and data["epg_normed"].values[0] is not None:
                epg = data["epg_normed"].values[0]
                self.fig_epgs[repl].update_layout(
                    yaxis_title="Signal intensity; (normalized wrt control peak)"
                )
            else:
                epg = data["epg_raw"].values[0]
            nucleotides = epg.nucleotides
            trace = epg.trace
            left_bound = data.left_bound.values[0]
            right_bound = data.right_bound.values[0]
            fail_reason = data.failed.values[0]
            label = f"t={str(timepoint) + 'h': <12}{'[FAILED] ' + fail_reason if fail_reason else ''}"
            # add epg
            self.fig_epgs[repl].add_trace(
                go.Scatter(
                    x=nucleotides,
                    y=trace,
                    name=label,
                    mode="lines",
                    line=dict(
                        width=3,
                        color=self.timepoint_to_color[timepoint].format(opacity=0.95),
                    ),
                ),
            )
            if left_bound is not None and fail_reason is None:
                if timepoint == 0:
                    tp0_left_bound = left_bound
                    tp0_right_bound = right_bound
                    # add target bounds
                    self.fig_epgs[repl].add_vline(
                        tp0_left_bound,
                        line=dict(
                            dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"
                        ),
                    )
                    self.fig_epgs[repl].add_vline(
                        tp0_right_bound,
                        line=dict(
                            dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"
                        ),
                    )
                # add area under the curve of target
                self.fig_epgs[repl].add_traces(
                    self._get_area_traces(
                        nucleotides=nucleotides,
                        trace=trace,
                        boundaries=(tp0_left_bound, tp0_right_bound),
                        color=self.timepoint_to_color[timepoint],
                        remove_background=self.remove_background,
                        prominence=data.prominence.values[0],
                    )
                )

            if self.control_peak_min is not None and timepoint == 0:
                # add control bounds
                self.fig_epgs[repl].add_vline(
                    self.control_peak_min,
                    line=dict(dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"),
                )
                self.fig_epgs[repl].add_vline(
                    self.control_peak_max,
                    line=dict(dash="dash", width=1, color="rgba(249, 249, 249, 0.5)"),
                )

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

    def combined_plotting(self) -> None:
        self.plot_decay_curve()
        self.plot_residuals()
        self.plot_rna_amounts()
        self.plot_sample_traces()

    @staticmethod
    def _get_area_traces(
        nucleotides: npt.NDArray[np.float64],
        trace: npt.NDArray[np.float64],
        boundaries: tuple[float, float],
        color: str,
        remove_background: bool,
        prominence: float,
    ) -> go.Scatter:
        mask = (nucleotides > boundaries[0]) & (nucleotides < boundaries[1])
        nucleotides = nucleotides[mask].copy()
        trace = trace[mask].copy()
        if remove_background:
            base = trace.max() - prominence
            background_line = [base] * len(trace)
        else:
            background_line = [0] * len(trace)
        y = background_line + trace[::-1].tolist()
        x = nucleotides.tolist() + nucleotides[::-1].tolist()
        return [
            go.Scatter(
                x=x,
                y=y,
                fill="toself",
                mode="lines",
                fillcolor=color.format(opacity=0.35),
                line_width=0,
                showlegend=False,
            ),
            go.Scatter(
                x=nucleotides,
                y=background_line,
                mode="lines",
                line=dict(color=color.format(opacity=0.75), dash="dash", width=3),
                showlegend=False,
            ),
        ]

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
