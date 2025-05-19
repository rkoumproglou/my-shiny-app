from shiny import App, ui, render, reactive
import pandas as pd
from collections import Counter
from scipy.stats import chisquare
import plotly.graph_objects as go

# Define Mendelian models
MENDELIAN_MODELS = {
    "3:1": [3, 1],
    "1:2:1": [1, 2, 1],
    "1:1:1:1": [1, 1, 1, 1],
    "9:3:3:1": [9, 3, 3, 1],
    "12:3:1": [12, 3, 1],
    "9:7": [9, 7],
    "15:1": [15, 1],
    "9:3:4": [9, 3, 4],
    "13:3": [13, 3],
    "9:6:1": [9, 6, 1]
}

# Interpretation text for each model
MODEL_INTERPRETATIONS = {
    "3:1": "The model 3:1 indicates **monogenic dominance**. In this case, a single gene has a dominant allele that masks the effect of a recessive one. This results in three offspring showing the dominant trait for every one showing the recessive trait.",
    "1:2:1": "The model 1:2:1 indicates **additive gene action**. In this case, one gene with two alleles produces three phenotypes where the heterozygote expresses an intermediate form.",
    "12:3:1": "The model 12:3:1 indicates **dominant epistasis**. In this case, a dominant allele at one gene locus masks the expression of another gene, reducing the number of phenotypes to three.",
    "9:7": "The model 9:7 indicates **duplicate recessive epistasis**. In this case, both genes must have at least one dominant allele to produce a specific phenotype, otherwise a common phenotype appears.",
    "15:1": "The model 15:1 indicates **duplicate dominant epistasis**. In this case, a dominant allele at either of two loci produces the same phenotype, resulting in a very rare unique phenotype.",
    "9:3:4": "The model 9:3:4 indicates **recessive epistasis**. In this case, one gene suppresses the expression of another gene only when in homozygous recessive form.",
    "13:3": "The model 13:3 indicates **dominant and recessive (inhibitory) epistasis**. In this case, the presence of one dominant allele inhibits the phenotypic expression of another gene.",
    "9:6:1": "The model 9:6:1 indicates **polymeric gene interaction**. In this case, two genes contribute equally to a new phenotype, while single gene expressions produce intermediate phenotypes.",
    "1:1:1:1": "The model 1:1:1:1 indicates **independent assortment of two genes**. In this case, each combination of alleles is equally likely, often seen in a dihybrid test cross.",
    "9:3:3:1": "The model 9:3:3:1 indicates **independent assortment of two genes with complete dominance**. In this case, the classic Mendelian dihybrid ratio arises from two genes segregating independently."
}

# Define UI
app_ui = ui.page_fluid(
    ui.h2("Mendelian Ratio Chi-square Tester"),
    ui.layout_columns(
        ui.card(
            ui.markdown("**Paste your phenotypic data (one value per line):**"),
            ui.input_text_area("phenodata", "", rows=15),
        ),
        #ui.hr(),
        ui.card(
            ui.output_ui("result_ui"),
        ),
        col_widths=(3,9)
    )
)

# Define server logic
def server(input, output, session):

    @reactive.calc
    def parsed_data():
        lines = input.phenodata().strip().splitlines()
        cleaned = [line.strip().lower() for line in lines]
        return Counter(cleaned)

    @reactive.calc
    def matched_models():
        count = len(parsed_data())
        return {
            name: ratio
            for name, ratio in MENDELIAN_MODELS.items()
            if len(ratio) == count
        }

    @reactive.calc
    def test_results():
        obs_counts = parsed_data()
        total = sum(obs_counts.values())
        sorted_categories = sorted(obs_counts.keys())

        results = []

        for name, ratio in matched_models().items():
            expected = sorted([total * r / sum(ratio) for r in ratio])
            observed = sorted([obs_counts.get(cat, 0) for cat in sorted_categories])

            if len(observed) != len(expected):
                continue
            chi2, p = chisquare(f_obs=observed, f_exp=expected)
            results.append((name, chi2, p, observed, expected))

        return sorted(results, key=lambda x: -x[2])

    @output
    @render.ui
    def result_ui():
        if not parsed_data():
            return ui.p("Paste data to begin.")

        if not matched_models():
            return ui.p("❌ No Mendelian model matches the number of phenotypic categories.")

        best = test_results()[0]
        name, chi2, p, obs, exp = best
        categories = sorted(parsed_data().keys())

        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=categories,
            y=obs,
            name="Observed",
            marker_color="lightgreen"
        ))
        fig.add_trace(go.Bar(
            x=categories,
            y=exp,
            name=f"Expected ({name})",
            marker_color="darkgreen"
        ))
        fig.update_layout(
            barmode="group",
            title="Observed vs Expected Counts",
            legend_title="Model Explanation",
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)"
        )

        return ui.div(
            ui.card(
                ui.div(
                    ui.h4("Best Fit Model", style="color: darkgreen;"),
                    ui.p(f"Model: **{name}**"),
                    ui.p(f"Observed counts: {obs}"),
                    ui.p(f"Expected counts: {[round(e, 2) for e in exp]}"),
                    ui.p(f"Chi-square statistic: {chi2:.4f}"),
                    ui.p(f"P-value: {p:.4f}"),
                ),
                style="box-shadow: 2px 2px 10px #ccc; padding: 1rem;"
            ),
            ui.card(
                ui.div(
                    ui.h4("Model Interpretation", style="color: darkgreen;"),
                    ui.markdown(MODEL_INTERPRETATIONS.get(name, "No interpretation available for this model."))
                ),
                style="box-shadow: 2px 2px 10px #ccc; padding: 1rem; margin-top: 1rem;"
            ),
            ui.output_plot("bar_plot"),
        )

    @output
    @render.plot
    def bar_plot():
        if not test_results():
            return None
        name, _, _, obs, exp = test_results()[0]
        categories = sorted(parsed_data().keys())
        df = pd.DataFrame({
            "Category": categories,
            "Observed": obs,
            "Expected": exp
        })
        return df.plot.bar(x="Category", color=["lightgreen", "darkgreen"], rot=0).get_figure()

# Create and run the app
app = App(app_ui, server)
