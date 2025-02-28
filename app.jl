using CSV
using DataFrames
using QuartetNetworkGoodnessFit
using PhyloNetworks
using Distributions
using GenieFramework
using PlotlyBase
using Stipple
@genietools

const species_tree = readTopology("(Smi165:18.12505298,(Age001:8.747492442,(Adi001:8.308140028,(((Asu001:3.171531094,Aga001:3.171531094):2.767496398,(Ape001:1.381629136,(Aza037:1.363263922,Ama006:1.363263922):0.01836521346):4.555985647):0.1963665132,Aru001:6.135394005):2.172746023):0.4393524143):9.377560535);")
const gene = joinpath("uploads", "gene_tree.rtf")
genetrees = readMultiTopology(gene)
q,t = countquartetsintrees(genetrees)
df = writeTableCF(q,t)
CSV.write("tableCF.csv", df)
const geneCF = readTableCF("tableCF.csv")
topologyMaxQPseudolik!(species_tree, geneCF);
df_long = fittedQuartetCF(geneCF, :long)

mean_Norm = 0
stddev = 1

x1_range = range(-3*stddev, 3*stddev, length=100)
y1 = pdf(Normal(mean_Norm, stddev), x1_range)

res = quarnetGoFtest!(species_tree, geneCF,true;quartetstat=:LRT,correction=:simulation,seed=100,nsim=100,verbose=false,keepfiles=false)

@handlers begin
    @out TFchoice = ["true", "false"]
    @out stat_choice = ["LRT", "Qlog", "pearson"]
    @out corr_choice = ["simulation", "none"]

    @in optimal_bl = ["true"]
    @in test_stat = ["LRT"]
    @in correct = ["simulation"]
    @in setseed = 100
    @in nsimulation = 100
    @in diagnose = ["false"]
    @in keep = ["false"]
    
    @out overall_p = res[1]
    @out trace_cf = [PlotlyBase.scatter(x=df_long[:, :obsCF],y=df_long[:, :expCF], mode = "markers")]
    @out layout_cf = PlotlyBase.Layout(title="CF Scatter Plot", 
    xaxis=attr(title="Observed CF"), yaxis=attr(title="Expected CF"))
    @out trace_p = [PlotlyBase.histogram(x=res[4], name="P Value Distribution")]
    @out layout_p = PlotlyBase.Layout(title="P-value Distribution",
        xaxis_title="P-value",
        yaxis_title="Count")
    @out Z_combined_trace = [PlotlyBase.histogram(x=res[6], name="Z-Value Distribution", histnorm="probability density", opacity=0.5),
    scatter(x=x1_range, y=y1, mode="lines", name="Standard Normal", line=attr(color="red", width=2)),
    scatter(x=range(-3*res[3], 3*res[3], length=100), y=pdf(Normal(mean_Norm, res[3]), range(-3*res[3], 3*res[3], length=100)), mode="lines", name="Fitted Normal", line=attr(color="blue", width=2))]
    @out Zhist_layout = PlotlyBase.Layout(
        title="Z-Value Distribution",
        xaxis_title="Z-value",
        yaxis_title="Density",
        barmode="overlay",
        legend=attr(
            orientation="h"
        )
    )

    @onchange optimal_bl,test_stat,correct,setseed,nsimulation,diagnose,keep begin
        optimal_bl_bool = optimal_bl[1] == "true"
        diagnose_bool = diagnose[1] == "true"
        keep_bool = keep[1] == "true"
        test_stat_value = Symbol(test_stat[1])
        correct_value = Symbol(correct[1])
        res1 = quarnetGoFtest!(species_tree, geneCF,optimal_bl_bool;quartetstat=test_stat_value,correction=correct_value,seed=setseed,nsim=nsimulation,verbose=diagnose_bool,keepfiles=keep_bool)
        overall_p = res1[1]
        sigma = res1[3]
        x2_range = range(-3*sigma, 3*sigma, length=100)
        y2 = pdf(Normal(mean_Norm, sigma), x2_range)
        trace_p = [PlotlyBase.histogram(x=res1[4], name="P Value Distribution")]
        Z_combined_trace = [PlotlyBase.histogram(x=res1[6], name="Z-Value Distribution", histnorm="probability density", opacity=0.5),
            scatter(x=x1_range, y=y1, mode="lines", name="Standard Normal", line=attr(color="red", width=2)),
            scatter(x=x2_range, y=y2, mode="lines", name="Fitted Normal", line=attr(color="blue", width=2))]
    end
end

# UI function
function ui()
    row([
    cell(class="st-col col-3", [
        h6("Set a seed"),
        slider(1:1:1000, :setseed),
        h6("Choose the number of simulations"),
        slider(1:1:1000, :nsimulation),
        Stipple.select(:optimal_bl, options = :TFchoice, label = "Optimal branch length"),
        Stipple.select(:test_stat, options = :stat_choice, label = "Choose a test statistics"),
        Stipple.select(:correct, options = :corr_choice, label = "Type for correction"),
        Stipple.select(:diagnose, options = :TFchoice, label = "Diagnose"),
        Stipple.select(:keep, options = :TFchoice, label = "Keep the analysis file"),
        p("The overall p-value is {{overall_p}}", class="st-module")
        ]),
    cell(class="st-col col-3", [
        plot(:trace_cf, layout=:layout_cf)
        ]),
    cell(class="st-col col-3", [
        plot(:trace_p, layout=:layout_p)
        ]),
    cell(class="st-col col-3", [
        plot(:Z_combined_trace, layout=:Zhist_layout)
        ])
    ])
end

@page("/",ui)