#nclude("Data.jl")
#include("OptMinConv.jl")
#include("OptMoreColsConv.jl")
#include("PreMoreColsConv.jl")


using PlotlyJS

## Plot of w_j over \ell, per competence k
x = 1:noOfDays
#res = obtainResults(20, 0, c[1], c[2])
#@elapsed obtainResults(20, 0, lambda)
w = r[1]
yvals1 = [sum(w[i][1,:]) for i in x]
yvals2 = [sum(w[i][2,:]) for i in x]
yvals3 = [sum(w[i][3,:]) for i in x]
yvals4 = [sum(w[i][4,:]) for i in x]
yvals5 = [sum(w[i][5,:]) for i in x]

plot([
    scatter(name="Competence 1", x=x, y=vec(yvals1), mode="markers+lines", line_width=1),
    scatter(name="Competence 2",x=x, y=vec(yvals2), mode="markers+lines", line_width=1),
    scatter(name="Competence 3",x=x, y=vec(yvals3), mode="markers+lines", line_width=1),
    scatter(name="Competence 4",x=x, y=vec(yvals4), mode="markers+lines", line_width=1),
    scatter(name="Competence 5",x=x, y=vec(yvals5), mode="markers+lines", line_width=1)
    ],
    Layout(title="Number of agents required per day", xaxis_title="Day", yaxis_title="Number of agents"))
  
##



## Bar chart over sum(a_hj*w_j) in comparison to sum(D_h)
a_hjl = r[2]
w_jkl = r[1]
#w_jkl = v
#D = ceil.(Int, r[4])
D=r[4]
val1 = zeros(1, 20)
val2 = zeros(1, 20)
val3 = zeros(1, 20)
val4 = zeros(1, 20)
val5 = zeros(1, 20)
D_h1 = zeros(1, 20)
D_h2 = zeros(1, 20)
D_h3 = zeros(1, 20)
D_h4 = zeros(1, 20)
D_h5 = zeros(1, 20)

for l in 1:size(a_hjl,1), h in 1:noOfTimeSteps, j in 1:size(a_hjl[l],2)
    val1[l] = val1[l] + a_hjl[l][h,j] * w_jkl[l][1,j]
end
for l in 1:size(a_hjl,1), h in 1:noOfTimeSteps, j in 1:size(a_hjl[l],2)
    val2[l] = val2[l] + a_hjl[l][h,j] * w_jkl[l][2,j]
end
for l in 1:size(a_hjl,1), h in 1:noOfTimeSteps, j in 1:size(a_hjl[l],2)
    val3[l] = val3[l] + a_hjl[l][h,j] * w_jkl[l][3,j]
end
for l in 1:size(a_hjl,1), h in 1:noOfTimeSteps, j in 1:size(a_hjl[l],2)
    val4[l] = val4[l] + a_hjl[l][h,j] * w_jkl[l][4,j]
end
for l in 1:size(a_hjl,1), h in 1:noOfTimeSteps, j in 1:size(a_hjl[l],2)
    val5[l] = val5[l] + a_hjl[l][h,j] * w_jkl[l][5,j]
end

for l in 1:noOfDays
    D_h1[l] = sum(D[l,:,1])
end
for l in 1:noOfDays
    D_h2[l] = sum(D[l,:,2])
end
for l in 1:noOfDays
    D_h3[l] = sum(D[l,:,3])
end
for l in 1:noOfDays
    D_h4[l] = sum(D[l,:,4])
end
for l in 1:noOfDays
    D_h5[l] = sum(D[l,:,5])
end

x = 1:noOfDays
y11 = val1
y12 = val2
y13 = val3
y14 = val4
y15 = val5
y21 = D_h1
y22 = D_h2
y23 = D_h3
y24 = D_h4
y25 = D_h5

c1 = plot([ # har för mig att färgen på staplarna var lite whack
    bar(name="Total availability", x=1:20, y=vec(y11)),
    bar(name="Total demand", x=1:20, y=vec(y21))
], 
Layout(title="Total availability/demand for competence 1", xaxis_title="Day", yaxis_title="Number of agents"))
relayout!(c1, barmode="group")

c2 = plot([
    bar(name="Total availability", x=1:20, y=vec(y12)),
    bar(name="Total demand", x=1:20, y=vec(y22))
], 
Layout(title="Total availability/demand for competence 2", xaxis_title="Day", yaxis_title="Number of agents"))
relayout!(c2, barmode="group")

c3 = plot([
    bar(name="Total availability", x=1:20, y=vec(y13)),
    bar(name="Total demand", x=1:20, y=vec(y23))
], 
Layout(title="Total availability/demand for competence 3", xaxis_title="Day", yaxis_title="Number of agents"))
relayout!(c3, barmode="group")

c4 = plot([
    bar(name="Total availability", x=1:20, y=vec(y14)),
    bar(name="Total demand", x=1:20, y=vec(y24))
], 
Layout(title="Total availability/demand for competence 4", xaxis_title="Day", yaxis_title="Number of agents"))
relayout!(c4, barmode="group")

c5 = plot([
    bar(name="Total availability", x=1:20, y=vec(y15)),
    bar(name="Total demand", x=1:20, y=vec(y25))
], 
Layout(title="Total availability/demand for competence 5", xaxis_title="Day", yaxis_title="Number of agents"))
relayout!(c5, barmode="group")

savefig(c1, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\PrebarAv071.pdf")
savefig(c2, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\barAv0952RMP.pdf")
savefig(c3, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\barAv0953RMP.pdf")
savefig(c4, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\barAv0954RMP.pdf")
savefig(c5, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\barAv0955RMP.pdf")
##

## Plot y_h/w_j over D_h
day = 1
#res = obtainResults(1, day, c[1], c[2])
#D = ceil.(Int, r[4])
D = r[4]
yu_lhk = r[3]
#w_kj = r[1][day]
#yu_lhk = yfrac[day]
w_kj = v[day]
a_hj = r[2][day]

xvals1 = D[day,:,1]
xvals2 = D[day,:,2]
xvals3 = D[day,:,3]
xvals4 = D[day,:,4]
xvals5 = D[day,:,5]

yvals1 = yu_lhk[day][:,1]
yvals2 = yu_lhk[day][:,2]
yvals3 = yu_lhk[day][:,3]
yvals4 = yu_lhk[day][:,4]
yvals5 = yu_lhk[day][:,5]

yvals11 = sum(w_kj[1,j] .* a_hj[:,j] for j in 1:size(w_kj[1,:],1))
yvals12 = sum(w_kj[2,j] .* a_hj[:,j] for j in 1:size(w_kj[2,:],1))
yvals13 = sum(w_kj[3,j] .* a_hj[:,j] for j in 1:size(w_kj[3,:],1))
yvals14 = sum(w_kj[4,j] .* a_hj[:,j] for j in 1:size(w_kj[4,:],1))
yvals15 = sum(w_kj[5,j] .* a_hj[:,j] for j in 1:size(w_kj[5,:],1))

yp = plot([
    bar(name="y_h",x=1:noOfTimeSteps, y=yvals1),
    bar(name="D_h",x=1:noOfTimeSteps, y=xvals1)],

    Layout(title="Availability/demand per time-step for competence 1, day 5", xaxis_title="h", yaxis_title="Number of agents")
)
relayout!(yp, barmode="group")

wp = plot([
    bar(name="Availability", x=1:noOfTimeSteps, y=yvals11),
    bar(name="Demand",x=1:noOfTimeSteps, y=xvals1)],

    Layout(title="Availability/demand per time step for competence 1, day 1", xaxis_title="Time step", yaxis_title="Number of agents")
)
relayout!(wp, barmode="group")

savefig(wp, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\barHRMP.pdf")
##

#show(yu_lhk[1][:,1]')
#w_kj[1,:]'

## Plot of c1*sum(y_h) over c2*sum(a_hj*w_j),   lambda is included in "Data.jl"
day = 10
cost = [1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]
xyvals = zeros(2, length(cost))
ywvals = zeros(2, length(cost))

for i in 1:length(cost)
    println("lambda = ", cost[i])
    res = obtainResults(20, 0, cost[i])
    yu_lhk = res[3]
    a_lhj = res[2]
    w_lkj = res[1]
    #xyvals[1,i] = (1-cost[i])*(sum(sum(sum((a_lhj[day][h,j])*w_lkj[day][k,j] for h in TIMESTEP_H) for j in 1:size(a_lhj[day],2) for k in COMPETENCE_K)) + 
    #sum(sum(((findlast(x -> x == 1, BitArray(a_lhj[day][:,j])) - findnext(BitArray(a_lhj[day][:,j]), 1)) >= 24 ? 2*(1-cost[i])*w_lkj[day][k,j] : (1-cost[i])*w_lkj[day][k,j]) for k in COMPETENCE_K) for j in 1:size(a_lhj[day],2)))
    #xyvals[2,i] = cost[i]*(sum(sum(yu_lhk[day][h,k] for h in TIMESTEP_H) for k in COMPETENCE_K))
    xyvals[1,i] = sum(sum(sum((a_lhj[day][h,j])*w_lkj[day][k,j] for h in TIMESTEP_H) for j in 1:size(a_lhj[day],2) for k in COMPETENCE_K))
    xyvals[2,i] = sum(sum(yu_lhk[day][h,k] for h in TIMESTEP_H) for k in COMPETENCE_K)
end

plot([
    scatter(name="lambda=1", x=[xyvals[1,1]], y=[xyvals[2,1]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.9", x=[xyvals[1,2]], y=[xyvals[2,2]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.8", x=[xyvals[1,3]], y=[xyvals[2,3]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.7", x=[xyvals[1,4]], y=[xyvals[2,4]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.6", x=[xyvals[1,5]], y=[xyvals[2,5]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.5", x=[xyvals[1,6]], y=[xyvals[2,6]], mode="markers+lines", line_width=1)],
    #scatter(name="lambda=0.4", x=[xyvals[1,7]], y=[xyvals[2,7]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.3", x=[xyvals[1,8]], y=[xyvals[2,8]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.2", x=[xyvals[1,9]], y=[xyvals[2,9]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.1", x=[xyvals[1,10]], y=[xyvals[2,10]], mode="markers+lines", line_width=1)],
    #scatter(name="lambda=0", x=[xyvals[1,11]], y=[xyvals[2,11]], mode="markers+lines", line_width=1)], 

    Layout(title="Pareto front", xaxis_title="(1-lambda)-term", yaxis_title="lambda-term")
)

plot([
    #scatter(name="lambda=1", x=[xyvals[1,1]], y=[xyvals[2,1]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.95", x=[xyvals[1,2]], y=[xyvals[2,2]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.9", x=[xyvals[1,3]], y=[xyvals[2,3]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.85", x=[xyvals[1,4]], y=[xyvals[2,4]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.8", x=[xyvals[1,5]], y=[xyvals[2,5]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.75", x=[xyvals[1,6]], y=[xyvals[2,6]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.7", x=[xyvals[1,7]], y=[xyvals[2,7]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.65", x=[xyvals[1,8]], y=[xyvals[2,8]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.6", x=[xyvals[1,9]], y=[xyvals[2,9]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.55", x=[xyvals[1,10]], y=[xyvals[2,10]], mode="markers+lines", line_width=1),
    scatter(name="lambda=0.5", x=[xyvals[1,11]], y=[xyvals[2,11]], mode="markers+lines", line_width=1)],
    #scatter(name="lambda=0.45", x=[xyvals[1,12]], y=[xyvals[2,12]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.4", x=[xyvals[1,13]], y=[xyvals[2,13]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.35", x=[xyvals[1,14]], y=[xyvals[2,14]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.3", x=[xyvals[1,15]], y=[xyvals[2,15]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.25", x=[xyvals[1,16]], y=[xyvals[2,16]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.2", x=[xyvals[1,17]], y=[xyvals[2,17]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.15", x=[xyvals[1,18]], y=[xyvals[2,18]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.1", x=[xyvals[1,19]], y=[xyvals[2,19]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0.05", x=[xyvals[1,20]], y=[xyvals[2,20]], mode="markers+lines", line_width=1),
    #scatter(name="lambda=0", x=[xyvals[1,21]], y=[xyvals[2,21]], mode="markers+lines", line_width=1)], 

    Layout(title="Pareto front for day 10", xaxis_title="Total availability", yaxis_title="Amount of understaffing")
)

a = plot([
    scatter(name="lambda=0.95", x=[xyvals[1,1]], y=[xyvals[2,1]], mode="markers"),
    scatter(name="lambda=0.9", x=[xyvals[1,2]], y=[xyvals[2,2]], mode="markers"),
    scatter(name="lambda=0.85", x=[xyvals[1,3]], y=[xyvals[2,3]], mode="markers"),
    scatter(name="lambda=0.8", x=[xyvals[1,4]], y=[xyvals[2,4]], mode="markers"),
    scatter(name="lambda=0.75", x=[xyvals[1,5]], y=[xyvals[2,5]], mode="markers"),
    scatter(name="lambda=0.7", x=[xyvals[1,6]], y=[xyvals[2,6]], mode="markers"),
    scatter(name="lambda=0.65", x=[xyvals[1,7]], y=[xyvals[2,7]], mode="markers"),
    scatter(name="lambda=0.6", x=[xyvals[1,8]], y=[xyvals[2,8]], mode="markers"),
    scatter(name="lambda=0.55", x=[xyvals[1,9]], y=[xyvals[2,9]], mode="markers"),
    scatter(name="lambda=0.5", x=[xyvals[1,10]], y=[xyvals[2,10]], mode="markers")],

    Layout(title="Pareto front", xaxis_title="(1-lambda)-term", yaxis_title="lambda-term")
)
#xyvals[1,i] = (1-cost[i])*(sum(sum(sum((a_lhj[1][h,j])*w_lkj[1][k,j] for h in TIMESTEP_H) for j in 1:size(a_lhj[1],2) for k in COMPETENCE_K)) + 
#    sum(sum(((findlast(x -> x == 1, BitArray(a_lhj[1][:,j])) - findnext(BitArray(a_lhj[1][:,j]), 1)) >= 24 ? 2*(1-cost[i])*w_lkj[1][k,j] : (1-cost[i])*w_lkj[1][k,j]) for k in COMPETENCE_K) for j in 1:size(a_lhj[1],2)))
 #   xyvals[2,i] = cost[i]*(sum(sum(yu_lhk[1][h,k] for h in TIMESTEP_H) for k in COMPETENCE_K))

b = plot([
    scatter(name="lambda=0.95-0.05", x=xyvals[1,:], y=xyvals[2,:], mode="markers+lines", line_width=1)],
    Layout(title="Pareto front", xaxis_title="(1-lambda)-term", yaxis_title="lambda-term")
    )


## bar plot of x_ijl summed over l and j to get total shifts taken by agent i during planning period
r = obtainResults(10, 0, 0.75)
finalSolution = MILPmodel(r[1], r[2], 10)
f1 = round.(Int, finalSolution[1])

x_i = []
for i in 1:size(AGENT_I, 1)
    x_i = push!(x_i, sum(sum(f1[i,j,l] for j in 1:size(r[2][l], 2)) for l in 1:20))
end

plot([
    bar(name="Total number of shifts", x=1:45, y=x_i)],
    Layout(title="Total number of shifts taken per agent during planning period", xaxis_title="Agent", yaxis_title="Total number of shifts")
)


# bar plot of sum_j w_jk
day = 1
#w_kj = r[1][day]
w_kj = v[day]
a_hj = r[2][day]
telia_w1 = TeliaDemand(1)[2][day,:]
telia_w2 = TeliaDemand(2)[2][day,:]
telia_w3 = TeliaDemand(3)[2][day,:]
telia_w4 = TeliaDemand(4)[2][day,:]
telia_w5 = TeliaDemand(5)[2][day,:]

w_l1 = sum(w_kj[1,j] .* a_hj[:,j] for j in 1:size(w_kj[1,:], 1))
w_l2 = sum(w_kj[2,j] .* a_hj[:,j] for j in 1:size(w_kj[2,:], 1))
w_l3 = sum(w_kj[3,j] .* a_hj[:,j] for j in 1:size(w_kj[3,:], 1))
w_l4 = sum(w_kj[4,j] .* a_hj[:,j] for j in 1:size(w_kj[4,:], 1))
w_l5 = sum(w_kj[5,j] .* a_hj[:,j] for j in 1:size(w_kj[5,:], 1))

forecast_1 = TeliaDemand(1)[1][day,:]
forecast_2 = TeliaDemand(2)[1][day,:]
forecast_3 = TeliaDemand(3)[1][day,:]
forecast_4 = TeliaDemand(4)[1][day,:]
forecast_5 = TeliaDemand(5)[1][day,:]

p1=plot([
    bar(name="Model", x=1:38, y=w_l1),
    bar(name="Telia", x=1:38, y=telia_w1),
    bar(name="Forecast calls", x=1:38, y=forecast_1)],
    Layout(title="Availability per time step, Model vs Telia, competence 1, day 1", xaxis_title="Time step", yaxis_title="Number of agents/calls")
)
p2=plot([
    bar(name="Model", x=1:38, y=w_l2),
    bar(name="Telia", x=1:38, y=telia_w2),
    bar(name="Forecast calls", x=1:38, y=forecast_2)],
    Layout(title="Availability per time step, Model vs Telia, competence 2, day 1", xaxis_title="Time step", yaxis_title="Number of agents/calls")
)

p3=plot([
    bar(name="Model", x=1:38, y=w_l3),
    bar(name="Telia", x=1:38, y=telia_w3),
    bar(name="Forecast calls", x=1:38, y=forecast_3)],
    Layout(title="Availability per time step, Model vs Telia, competence 3, day 1", xaxis_title="Time step", yaxis_title="Number of agents/calls")
)

p4=plot([
    bar(name="Model", x=1:38, y=w_l4),
    bar(name="Telia", x=1:38, y=telia_w4),
    bar(name="Forecast calls", x=1:38, y=forecast_4)],
    Layout(title="Availability per time step, Model vs Telia, competence 4, day 1", xaxis_title="Time step", yaxis_title="Number of agents/calls")
)

p5=plot([
    bar(name="Model", x=1:38, y=w_l5),
    bar(name="Telia", x=1:38, y=telia_w5),
    bar(name="Forecast calls", x=1:38, y=forecast_5)],
    Layout(title="Availability per time step, Model vs Telia, competence 5, day 1", xaxis_title="Time step", yaxis_title="Number of agents/calls")
)

savefig(p1, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\MvT0951.pdf")
savefig(p2, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\MvT0952.pdf")
savefig(p3, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\MvT0953.pdf")
savefig(p4, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\MvT0954.pdf")
savefig(p5, "C:\\Users\\Asus\\OneDrive\\Skrivbord\\JuliaPlots\\MvT0955.pdf")