
using Plots
using Distributions
using KernelDensity
""
function plotHist(samples,name::String, title::String)
    histogram(samples,normalize=:pdf,xlabel=name,title=title)
end

function plotPDF(samples,name::String, title::String)
    lo   = minimum(samples)
    high = maximum(samples)
    y=(high-lo)/100
    if y==0
        z=0.01
        x=collect(lo-3*z:z:high+3*z)
        plot!([lo],seriestype="vline",xlabel=name,title=title, label="")
    else
        x=collect(lo-3*y:y:high+3*y)
        plot(x,pdf(kde(samples),x),xlabel=name,title=title, label="")
    end
    
end

function plotHist_volt(result,name::String, mop, samplesize::Int=1000; busno= [], cdf=false, pdf=false) #the variable of interest
   no_coeff= length(result["solution"]["nw"])
   if length(busno) == 0
   bus_len=length(result["solution"]["nw"]["1"]["bus"])
   busno=[parse(Int64, r) for r in keys(result["solution"]["nw"]["1"]["bus"])]
   vx_coeff= [([(result["solution"]["nw"]["$i"]["bus"]["$j"][name]) for i in 1:no_coeff]) for j in keys(result["solution"]["nw"]["1"]["bus"])]
   else
    
    bus_len=length(busno)
   vx_coeff= [([(result["solution"]["nw"]["$i"]["bus"]["$j"][name]) for i in 1:no_coeff]) for j in busno]
    
end
   p_arr1 = Plots.Plot{Plots.GRBackend}[] 
   p_arr2 = Plots.Plot{Plots.GRBackend}[] 
   p_arr3 =  Plots.Plot{Plots.GRBackend}[]
 for i in 1:length(busno)
    z=_PCE.samplePCE(samplesize,vx_coeff[i],mop) #samples
    push!(p_arr1, plotHist(z,name, "bus_$(busno[i])"))
    push!(p_arr2, plot(sort(z),(1:samplesize)./samplesize,
    xlabel = name, ylabel = "Probability", 
    title = "Bus_$(busno[i]) $name CDF", label = ""))
    push!(p_arr3, plotPDF(z,name, "PDF bus_$(busno[i]) $name"))
 end
    display(plot(p_arr1...))

    if pdf==true
        display(plot(p_arr3...))
    end

    if cdf==true
        display(plot(p_arr2...))
    end
end

function plotHist_gen(result,name::String, mop, samplesize::Int=1000; cdf=false, pdf=false) #the variable of interest
    no_coeff= length(result["solution"]["nw"])
    gen_len=length(result["solution"]["nw"]["1"]["gen"])
    vx_coeff= [([(result["solution"]["nw"]["$i"]["gen"]["$j"][name]) for i in 1:no_coeff]) for j in keys(result["solution"]["nw"]["1"]["gen"])]
    p_arr1 = Plots.Plot{Plots.GRBackend}[] 
    p_arr2 = Plots.Plot{Plots.GRBackend}[] 
    p_arr3 =  Plots.Plot{Plots.GRBackend}[]
    for i in 1:gen_len
        z=_PCE.samplePCE(samplesize,vx_coeff[i],mop) #samples
        push!(p_arr1, plotHist(z,name, "gen_$i"))
        push!(p_arr2, plot(sort(z),(1:samplesize)./samplesize,
        xlabel = name, ylabel = "Probability", 
        title = "Load_$i $name CDF", label = ""))
        push!(p_arr3, plotPDF(z,name, "PDF gen_$i $name"))
    end
    display(plot(p_arr1...))
     if cdf==true
       
        display(plot(p_arr2...))
     end
     if pdf==true
        display(plot(p_arr3...))
    end
 end


 function plotHist_load(result,name::String, mop, samplesize::Int=1000; cdf=false, pdf=false) #the variable of interest
    no_coeff= length(result["solution"]["nw"])
    load_len=length(result["solution"]["nw"]["1"]["load"])
    vx_coeff= [([(result["solution"]["nw"]["$i"]["load"]["$j"][name]) for i in 1:no_coeff]) for j in keys(result["solution"]["nw"]["1"]["load"])]
    p_arr1 = Plots.Plot{Plots.GRBackend}[]
    p_arr2 = Plots.Plot{Plots.GRBackend}[] 
    p_arr3 =  Plots.Plot{Plots.GRBackend}[]
    for i in 1:load_len
        z=_PCE.samplePCE(samplesize,vx_coeff[i],mop) #samples
        push!(p_arr1, plotHist(z,name, "Load_$i $name"))
        push!(p_arr2, plot(sort(z),(1:samplesize)./samplesize,
        xlabel = name, ylabel = "Probability", 
        title = "Load_$i $name CDF", label = ""))
        push!(p_arr3, plotPDF(z,name, "PDF load_$i $name"))
     end
     display(plot(p_arr1...))
     if cdf==true
       
        display(plot(p_arr2...))
     end
     if pdf==true
        display(plot(p_arr3...))
    end
 end

 function plotHist_branch(result,name::String, mop, samplesize::Int=1000; cdf=false, pdf=false) #the variable of interest
    no_coeff= length(result["solution"]["nw"])
    branch_len=length(result["solution"]["nw"]["1"]["branch"])
    vx_coeff= [([(result["solution"]["nw"]["$i"]["branch"]["$j"][name]) for i in 1:no_coeff]) for j in keys(result["solution"]["nw"]["1"]["branch"])]
    p_arr1 = Plots.Plot{Plots.GRBackend}[]
    p_arr2 = Plots.Plot{Plots.GRBackend}[]
    p_arr3 =  Plots.Plot{Plots.GRBackend}[]
    for i in 1:branch_len
        z=_PCE.samplePCE(samplesize,vx_coeff[i],mop) #samples
        push!(p_arr1, plotHist(z,name, "branch_$i $name"))
        push!(p_arr2, plot(sort(z),(1:samplesize)./samplesize,
        xlabel = name, ylabel = "Probability", 
        title = "Load_$i $name CDF", label = ""))
        push!(p_arr3, plotPDF(z,name, "PDF branch_$i $name")) 
     end
     
     display(plot(p_arr1...))
     if cdf==true 
        display(plot(p_arr2...))
     end
     if pdf==true
        display(plot(p_arr3...))
    end
 end