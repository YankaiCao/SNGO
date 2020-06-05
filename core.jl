const n_rel = 100
const lambda = 0.25
const alpha = 0.1
const mu_score = 0.15
const lambda_consecutive = 4
const sigma_violation = 1e-8
const con_tol = 1e-6
const mingap = 1e-2
const LP_improve_tol = min(1e-2, mingap)
const default_upper_bound_value = 1e8
const default_lower_bound_value = -1e8

const machine_error = 1e-10
const small_bound_improve = 1e-4
const large_bound_improve = 1e-1
const probing_improve = 1e-4
const nR = 0
const probing = true
const max_score_min = 1e-3
const OBBT = true
const debug = true
#const hasBin = false


function multi_start!(P, UB, n_trial = 3)
      updated = false
      lb = P.colLower
      ub = P.colUpper
      local_x = copy(P.colVal)
      P.solver = IpoptSolver(max_cpu_time = 600.0, print_level = 0)
      ncols = length(lb)
      mlb = -1e3*ones(ncols)
      mub =  1e3*ones(ncols)
      for i in 1:ncols
          if lb[i] != -Inf
             mlb[i] = lb[i]
          end
          if ub[i] != Inf
             mub[i] = ub[i]
	       end
      end
      for local_trial = 1:n_trial
	  if local_trial !=  1
	      percent = rand(ncols)
              initial  = mlb + (mub-mlb).*percent
              P.colVal = initial
	  end
          prime_status_trial = solve(P)
	  println("prime_status_trial   ",prime_status_trial, "    ", getobjectivevalue(P))
          if prime_status_trial == :Optimal
              local_obj_trial = getobjectivevalue(P)
              if local_obj_trial - UB < - mingap/10
                  UB = local_obj_trial
                  updated = true 
                  local_x = P.colVal
              end
          end
      end
      println("updated   ", updated, "   ", UB)
      local_status = :UnOptimal
      if updated
      	  local_status = :Optimal
      end
      P.colVal = copy(local_x)
      P.objVal = UB
      return (UB, local_status)
end




type StoSol
    firstSol::Vector{Float64}
    firstobjVal
    secondSols
    secondobjVals::Vector{Float64}
end
StoSol(nscen) = StoSol(Float64[], 0, Array{Array{Float64}}(nscen), Array{Float64}(nscen))


#=
function getStoSol(P::JuMP.Model)
    sol = StoSol()
    sol.firstSol = P.colVal
    sol.firstobjVal = P.objVal
    sol.secondSols = []
    sol.secondobjVals = []
    for (idx, scenario) in enumerate(PlasmoOld.getchildren(P))
        push!(sol.secondSols, scenario.colVal)
        push!(sol.secondobjVals, scenario.objVal)
    end
    return sol
end
=#

type Node
    xl::Vector{Float64}
    xu::Vector{Float64}
    bVarId::Int
    level::Int
    LB::Float64
    direction::Int    # -1 if it is left child, 1 if it is a right child, 0 if it is already used to updated the variable section, or if it is the root
    parWSSol
    parRSol
    parmaxScore::Float64
    colCat::Vector{Symbol}
end
Node() = Node(Float64[], Float64[], -1, -1, -1e10, 0, nothing, nothing, 0, Symbol[])

type PreprocessResult
     branchVarsId
     EqVconstr
     multiVariable_list
     multiVariable_convex
     multiVariable_aBB
     qbVarsId
     expVariable_list
     logVariable_list
     powerVariable_list
     monomialVariable_list 
end
PreprocessResult() = PreprocessResult(Int[], LinearConstraint[], [], [], [], [], [], [], [], [])



#lvar op exp(b * nlvar)
type ExpVariable
     lvarId::Int
     nlvarId::Int	    
     op     
     b::Float64
     cid::Vector{Int64}
end


#lvar op log(b * nlvar)
type LogVariable
     lvarId::Int
     nlvarId::Int
     op
     b::Float64
     cid::Vector{Int64}
end


# lvar op (b*nlvar)^d
type PowerVariable
     lvarId::Int
     nlvarId::Int
     op
     b::Float64
     d::Float64
     cid::Vector{Int64}
end

# lvar op (d)^(b*nlvar)
type MonomialVariable
     lvarId::Int
     nlvarId::Int
     op
     b::Float64
     d::Float64
     cid::Vector{Int64}
end



type MultiVariable
     terms::QuadExpr
     Q::Array{Float64,2}
     pd::Int   # =1 if convex, = -1 if concave, =0 if else
     qVarsId::Vector{Int}
     bilinearVarsId::Vector{Int}
     alpha::Vector{Float64}
end
MultiVariable()  = MultiVariable(zero(QuadExpr), zeros(Float64, 0, 0), 0, Int[], Int[], Float64[])


type MultiVariableCon
     mvs::Vector{MultiVariable}   #list of mvs
     aff::AffExpr    #remaining affine item
end

function Base.copy(mv::MultiVariable, new_model::Model)
    newmv = MultiVariable()
    newmv.terms = copy(mv.terms, new_model)
    newmv.Q = copy(mv.Q)
    newmv.pd = mv.pd
    newmv.qVarsId = copy(mv.qVarsId)
    newmv.bilinearVarsId = copy(mv.bilinearVarsId)
    newmv.alpha = copy(mv.alpha)
    newmv
end
function Base.copy(mvc::MultiVariableCon, new_model::Model)
    MultiVariableCon(map(x->copy(x, new_model), mvc.mvs), copy(mvc.aff, new_model))
end


type VarSelector
    varId::Vector{Int}
    score::Vector{Float64}
    pcost_right::Vector{Float64}
    pcost_left::Vector{Float64}
    n_right::Vector{Int}
    n_left::Vector{Int}
    tried::Vector{Bool}   # in the current node
end
VarSelector() = VarSelector(Int[], Float64[], Float64[], Float64[], Int[], Int[], Bool[])
VarSelector(bVarsId) = VarSelector(copy(bVarsId), ones(length(bVarsId)), zeros(length(bVarsId)), zeros(length(bVarsId)),zeros(Int, length(bVarsId)),zeros(Int, length(bVarsId)), falses(length(bVarsId)))
function sortVarSelector(b)
    perm = sortperm(b.score, rev=true)
    b.varId = b.varId[perm]
    b.score = b.score[perm]
    b.pcost_right = b.pcost_right[perm]
    b.pcost_left = b.pcost_left[perm]
    b.n_right = b.n_right[perm]
    b.n_left = b.n_left[perm]
end
function compute_score(left, right)
    score = (1-mu_score)*min(left, right)+mu_score*max(left, right)
end


function updateScore!(vs, Pex, P, Rsol, WSfirst, WS_status, relaxed_status)
    for i in 1:length(vs.varId)
        varId = vs.varId[i]
        bVarl = Pex.colLower[varId]
        bVaru = Pex.colUpper[varId]
        bValue = computeBvalue(Pex, P, varId, Rsol, WSfirst, WS_status, relaxed_status)
        nlinfeasibility_left = bValue - bVarl
        improve_left = vs.pcost_left[i]*nlinfeasibility_left
	nlinfeasibility_right = bVaru - bValue
        improve_right = vs.pcost_right[i]*nlinfeasibility_right
        vs.score[i] = compute_score(improve_left, improve_right)
        vs.tried[i] = false
    end
    sortVarSelector(vs)
end

function updateVarSelector(vs, node, P, node_LB)
    if node.direction != 0
        bvarId = node.bVarId
        i = findin(vs.varId, bvarId)[1]
        bVarl = P.colLower[bvarId]
        bVaru = P.colUpper[bvarId]
        improvement = node_LB - node.LB
        nlinfeasibility = bVaru - bVarl
        if node.direction == -1
            vs.pcost_left[i] = (vs.pcost_left[i]*vs.n_left[i] + improvement/nlinfeasibility)/(vs.n_left[i] + 1)
            vs.n_left[i] += 1
        elseif node.direction == 1
            vs.pcost_right[i] = (vs.pcost_right[i]*vs.n_right[i] + improvement/nlinfeasibility)/(vs.n_right[i] + 1)
            vs.n_right[i] += 1
        end
    end
end

#find the node with the lowest lower bound
function getGlobalLowerBound(nodeList, FLB, UB)
    LB = 1e10
    nodeid = 1
    for (idx,n) in enumerate(nodeList)
        if n.LB < LB
            LB = n.LB
            nodeid = idx
        end
    end
    #LB = min(LB, FLB)
    LB = min(LB, UB)
    return LB, nodeid
end

function getRobjective(R, relaxed_status, defaultLB)
    node_LB = defaultLB
    if relaxed_status == :Optimal
        node_LB = getobjectivevalue(R)
    elseif relaxed_status == :Infeasible
        #if relaxed problem is infeasible, delete
        node_LB = 2e10
    end
    return node_LB
end



function updateStoFirstBounds!(P)
    n = P.numCols	 
    for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
        firstVarsId = scenario.ext[:firstVarsId]
        for i in 1:n
	    secondId = firstVarsId[i] # the Id of second stage variable corresponds to first stage Id
            if secondId > 0			       		         
	        if hasBin
                    if scenario.colCat[secondId] == :Bin
                        if scenario.colLower[secondId] == 0.0 && scenario.colUpper[secondId] == 1.0
                            if P.colCat[i] == :Cont 
                                P.colCat[i] = :Bin
				if P.colLower[i] < 0.0
				    P.colLower[i] = 0.0
				end
				if P.colUpper[i] > 1.0
				    P.colUpper[i] == 1.0
				end    
                            end			    
                        elseif scenario.colLower[secondId] >= machine_error && scenario.colUpper[secondId] == 1.0
			    fixVar(Variable(scenario, secondId), 1.0)			
			    fixVar(Variable(P, i), 1.0)
			elseif scenario.colLower[secondId] == 0.0  && 1.0 - scenario.colUpper[secondId] >= machine_error			
			    fixVar(Variable(scenario, secondId), 0.0)
			    fixVar(Variable(P, i), 0.0)			    
                        end
                    elseif scenario.colCat[secondId] == :Fixed
		    	fixVar(Variable(P, i), scenario.colLower[secondId])   
		    end
                end
		
                if scenario.colLower[secondId] >  P.colLower[i]
                    P.colLower[i] = scenario.colLower[secondId]
                end
                if scenario.colUpper[secondId] <  P.colUpper[i]
                    P.colUpper[i] = scenario.colUpper[secondId]
                end
            end
        end
    end
    if hasBin
        updateFirstBounds!(P, P.colLower, P.colUpper, P.colCat)
    else
        updateFirstBounds!(P, P.colLower, P.colUpper)
    end
end

# update bounds and cat information in the first stage variables of the master-scenario form
function updateFirstBounds!(P, xl::Array{Float64,1}, xu::Array{Float64,1}, cat=nothing)
        P.colLower = copy(xl)
        P.colUpper = copy(xu)
	if cat != nothing
	    P.colCat = copy(cat)
	end

	n = P.numCols
        if hasBin
	    for i in 1:n
	        if P.colCat[i] == :Bin
                    if P.colLower[i] >= machine_error && P.colUpper[i] == 1.0
                        fixVar(Variable(P, i), 1.0)
                    elseif P.colLower[i] == 0.0  && 1.0 - P.colUpper[i] >= machine_error
                        fixVar(Variable(P, i), 0.0) 	       
                    end
		end	
	    end
	end
	#=
	for i = 1:n
	    if xl[i] == xu[i]  
	       	fixVar(Variable(P, i), xl[i])     
	    end   
	end
	=#
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))	
	    for i = 1:n
	    	secondId = scenario.ext[:firstVarsId][i] # the Id of second stage variable corresponds to first stage Id
	    	if secondId != -1
            	    scenario.colLower[secondId] = copy(xl[i])
           	    scenario.colUpper[secondId] = copy(xu[i])
		    if cat != nothing
 		        scenario.colCat[secondId] = P.colCat[i]
		    end
                end
	    end
        end
end

# update bounds of scenarios from master
function updateFirstBounds!(P, xl::Float64, xu::Float64, varId::Int64, cat=nothing)
        P.colLower[varId] = copy(xl)
        P.colUpper[varId] = copy(xu)
	if cat != nothing
	    P.colCat[varId] = cat
	end    
        if hasBin
            if P.colCat[varId] == :Bin
                if P.colLower[varId] >= machine_error && P.colUpper[varId] == 1.0
                    fixVar(Variable(P, varId), 1.0)
                elseif P.colLower[varId] == 0.0  && 1.0 - P.colUpper[varId] >= machine_error
                    fixVar(Variable(P, varId), 0.0)
                end
            end
        end
	#=
        if xl[varId] ==	xu[varId]
	   	fixVar(Variable(P, varId), xl[varId])     
	end
	=# 
	for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
	    secondId = scenario.ext[:firstVarsId][varId]
	    if secondId != -1
                scenario.colLower[secondId] = copy(xl)
                scenario.colUpper[secondId] = copy(xu)
		if cat != nothing
                    scenario.colCat[secondId] = P.colCat[varId]
		end
	    end   
        end
end

# update bounds information of extensive form
function updateExtensiveFirstBounds!(Pex, P, xl::Array{Float64,1}, xu::Array{Float64,1}, cat=nothing)
        n = P.numCols
        Pex.colLower[1:n] = copy(xl)
        Pex.colUpper[1:n] = copy(xu)
        if cat != nothing
            Pex.colCat[1:n] = copy(cat)
        end
        if hasBin
            for i in 1:n
                if Pex.colCat[i] == :Bin
                    if Pex.colLower[i] >= machine_error && Pex.colUpper[i] == 1.0
                        fixVar(Variable(Pex, i), 1.0)
                    elseif Pex.colLower[i] == 0.0  && 1.0 - Pex.colUpper[i] >= machine_error
                        fixVar(Variable(Pex, i), 0.0)
                    end
                end
            end
        end
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            for j = 1:length(scenario.ext[:firstVarsId])
	    	secondId = scenario.ext[:firstVarsId][j]
                if secondId != -1
            	    Pex.colLower[v_map[secondId]] = copy(Pex.colLower[j])
            	    Pex.colUpper[v_map[secondId]] = copy(Pex.colUpper[j])
		    if cat != nothing
		        Pex.colCat[v_map[secondId]] = Pex.colCat[j]
		    end
		end
	    end	    
        end
end

function updateExtensiveFirstBounds!(Pex, P, xl::Float64, xu::Float64, varId::Int64, cat=nothing)
        n = P.numCols
        Pex.colLower[varId] = copy(xl)
        Pex.colUpper[varId] = copy(xu)
        if cat != nothing
            Pex.colCat[varId] = copy(cat)
        end

        if hasBin
            if Pex.colCat[varId] == :Bin
                if Pex.colLower[varId] >= machine_error && Pex.colUpper[varId] == 1.0
                    fixVar(Variable(Pex, varId), 1.0)
                elseif Pex.colLower[varId] == 0.0  && 1.0 - Pex.colUpper[varId] >= machine_error
                    fixVar(Variable(Pex, varId), 0.0)
                end
            end
        end

        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
	    secondId = scenario.ext[:firstVarsId][varId]
	    if secondId != -1
                Pex.colLower[v_map[secondId]] = copy(Pex.colLower[varId])
            	Pex.colUpper[v_map[secondId]] = copy(Pex.colUpper[varId])
		if cat != nothing
		    Pex.colCat[v_map[secondId]] = Pex.colCat[varId]
		end
	    end	
        end
end

function updateExtensiveFirstBounds!(PexLower, PexUpper, P, Pex, xl, xu, varId, PexCat = nothing, cat = nothing)
        n = P.numCols
        PexLower[varId] = copy(xl)
        PexUpper[varId] = copy(xu)

        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
	    secondId = scenario.ext[:firstVarsId][varId]
	    if secondId != -1
                PexLower[v_map[secondId]] = copy(xl)
                PexUpper[v_map[secondId]] = copy(xu)
            end
        end
	if PexCat != nothing
            PexCat[varId] = cat
            for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            	v_map = Pex.ext[:v_map][idx]
		secondId = scenario.ext[:firstVarsId][varId]
            	if secondId != -1
                    PexCat[v_map[secondId]] = cat
                end
 	    end
        end
end

#=
function updateExtensiveFirstCat!(PexCat, Pex, cat, varId)
        PexCat[varId] = cat
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
	    secondId = scenario.ext[:firstVarsId][varId]
            if secondId != -1
                PexCat[v_map[secondId]] = cat
            end
        end
end
=#

## update bounds info from extensive form to master-children form
function updateExtensiveBoundsFromSto!(P, Pex)
        nfirst = P.numCols
        Pex.colLower[1:nfirst] = copy(P.colLower)
       	Pex.colUpper[1:nfirst] = copy(P.colUpper)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            Pex.colLower[v_map] = copy(scenario.colLower)
            Pex.colUpper[v_map] = copy(scenario.colUpper)
        end	
	if hasBin
	    Pex.colCat[1:nfirst] = copy(P.colCat) 
            for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            	v_map = Pex.ext[:v_map][idx]
            	Pex.colCat[v_map] = copy(scenario.colCat)
            end	    
	end
	
end

##bounds from extensive form to master-children form
function updateStoBoundsFromExtensive!(Pex, P)
	#updateStoFirstBounds!(P) 
        nfirst = P.numCols
        P.colLower = copy(Pex.colLower[1:nfirst])
        P.colUpper = copy(Pex.colUpper[1:nfirst])
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            scenario.colLower = copy(Pex.colLower[v_map])
            scenario.colUpper = copy(Pex.colUpper[v_map])
        end
	if hasBin
	    P.colCat = copy(Pex.colCat[1:nfirst])
            for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            	v_map = Pex.ext[:v_map][idx]
            	scenario.colCat = copy(Pex.colCat[v_map])
	    end	
        end
end


function updateStoSolFromSto!(P, dest)
        nfirst = P.numCols
        dest.colVal = copy(P.colVal[1:nfirst])
	dscenarios = PlasmoOld.getchildren(dest)

        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            dscenarios[idx].colVal = copy(scenario.colVal) 
            if in("objective_value", scenario.colNames)
	           dscenarios[idx].objVal = getvalue(scenario[:objective_value])
            end
        end
end


function updateStoSolFromExtensive!(Pex, P)
        nfirst = P.numCols
        P.colVal = copy(Pex.colVal[1:nfirst])
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            scenario.colVal = copy(Pex.colVal[v_map])
            if in("objective_value", scenario.colNames)   
               scenario.objVal = getvalue(scenario[:objective_value])
            end
        end
end

function SelectVarMaxRange(bVarsId, P::JuMP.Model, weights = ones(length(bVarsId)))
    bVarId = bVarsId[1]
    maxrange = 1e-8
    for v in bVarsId
        lb = P.colLower[v]
        if (lb == -Inf)
            lb = -1e8
        end
        ub = P.colUpper[v]
        if (ub == Inf)
            ub = 1e8
        end
        range = (ub - lb)*weights[v]
        if range > maxrange
            bVarId = v
            maxrange = range
        end
    end
    return bVarId
end

function updateStoBoundsFromSto!(m, dest)	
        dest.colLower = copy(m.colLower)
        dest.colUpper = copy(m.colUpper)
	dest.colCat = copy(m.colCat)
	mchildrens = PlasmoOld.getchildren(m)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(dest))
            scenario.colLower = copy(mchildrens[idx].colLower)
            scenario.colUpper = copy(mchildrens[idx].colUpper)
	    scenario.colCat = copy(mchildrens[idx].colCat)
        end
end

function inBounds(m::JuMP.Model, sol)
        sum(m.colLower.<=sol.<=m.colUpper) == length(sol)
end


function eval_g(c::AffExpr, x)
    sum = c.constant
    for i in 1:length(c.vars)
        sum += c.coeffs[i] * x[c.vars[i].col]
    end
    return sum
end

function eval_g(c::QuadExpr, x)
    sum = eval_g(c.aff, x)
    for i in 1:length(c.qcoeffs)
        sum += c.qcoeffs[i] * x[c.qvars1[i].col] * x[c.qvars2[i].col]
    end
    return sum
end

function eval_g(m::JuMP.Model, x)
    nl = length(m.linconstr)
    nq = length(m.quadconstr)
    n = nl + nq
    g = zeros(n, 1)
    for i in 1:nl
        g[i] = eval_g(m.linconstr[i].terms, x)
    end
    for i in 1:nq
        g[i+nl] = eval_g(m.quadconstr[i].terms, x)
    end
    return g
end


function computeBvalue(Pex, P, bVarId, Rsol, WSfirst, WS_status, relaxed_status)
    bVarl = Pex.colLower[bVarId]
    bVaru = Pex.colUpper[bVarId]
    mid = (bVarl+bVaru)/2
    bValue = mid
    if WS_status == :Optimal
        bValue = WSfirst[bVarId]
    elseif relaxed_status == :Optimal
        bValue = Rsol[bVarId]
    elseif UB_status == :Optimal
        bValue = P.colVal[bVarId]
    end
    if bValue >= bVaru || bValue <= bVarl
       bValue = mid
    end
    bValue = lambda*mid + (1-lambda)*bValue
    bValue = min(max(bValue, bVarl + alpha*(bVaru-bVarl)), bVaru - alpha*(bVaru-bVarl))
    return bValue
end


function projection!(sol, xl, xu)
    for i = 1:length(sol)    	 
    	if sol[i] <= xl[i]
	   sol[i] = xl[i]
	end   
	if sol[i] >= xu[i]
	   sol[i] = xu[i]
	end
    end
end


function factorable!(P)
    m = Model()
    m.solver = P.solver  
    # Variables
    m.numCols = P.numCols
    m.colNames = P.colNames[:]
    m.colNamesIJulia = P.colNamesIJulia[:]
    m.colLower = P.colLower[:]
    m.colUpper = P.colUpper[:]
    m.colCat = P.colCat[:]
    m.colVal = P.colVal[:]

    # Constraints
    m.linconstr  = map(c->copy(c, m), P.linconstr)
    m.quadconstr = map(c->copy(c, m), P.quadconstr)

    # Objective
    m.obj = copy(P.obj, m)
    m.objSense = P.objSense


    if !isempty(P.ext)
        m.ext = similar(P.ext)
        for (key, val) in P.ext
            m.ext[key] = try
                copy(P.ext[key])
            catch
                continue;  #error("Error copying extension dictionary. Is `copy` defined for all your user types?")
            end
        end
    end
 
    if P.nlpdata != nothing
        d = JuMP.NLPEvaluator(P)         #Get the NLP evaluator object.  Initialize the expression graph
        MathProgBase.initialize(d,[:ExprGraph])
        num_cons = MathProgBase.numconstr(P)
        for i = (1+length(m.linconstr)+length(m.quadconstr)):num_cons
            expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
            _modifycon!(m, expr)                  #splice the variables from v_map into the expression
        end
    end
    #=
    R2 = copyModel(m)
    d = JuMP.NLPEvaluator(R2)
    MathProgBase.initialize(d,[:ExprGraph])
    num_cons = MathProgBase.numconstr(R2)
    for i = 1:num_cons
            expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
            println(expr)
    end
    =#
    return m
end

function bounded(a)
    if abs(a)<1e8	 
        return true
    end
    return false
end



function positiveEven(a)
    if (a>0) && (a%2 == 0)
        return true
    end 
    return false	 
end

function negativeEven(a)
    if (a<0) && (a%2 == 0)
        return true
    end
    return false
end

function positiveOdd(a)
    if (a>0) && (a%2 == 1)
    return true
    end
    return false	
end

function negativeOdd(a)
    if (a<0) && (a%2 == -1)
        return true
    end
    return false
end

function Odd(a)
    if (a%2==-1) || (a%2 == 1)
        return true
    end
    return false
end

function positiveFrac(a)
    if (a>0) && (a%2 != 1) && (a%2 != 0) && (a%2 != -1)
    return true
    end
    return false	
end

function negativeFrac(a)
    if (a<0) && (a%2 != 1) && (a%2 != 0)&& (a%2 != -1)
        return true
    end
    return false
end



function scaleAff(a::AffExpr, c::Float64)
    a.coeffs = c*a.coeffs
    a.constant = c*a.constant 
    return a
end

function addAff(a::AffExpr, b::AffExpr)
    a.coeffs = [a.coeffs; b.coeffs]	 
    a.vars = [a.vars; b.vars]
    a.constant = a.constant + b.constant
    return a
end

function minusAff(a::AffExpr, b::AffExpr)
    a.coeffs = [a.coeffs; -b.coeffs]
    a.vars = [a.vars; b.vars]
    a.constant = a.constant - b.constant
    return a
end


function multAff(m, aff1, aff2, c)
     quad = QuadExpr()	 
     n1 = length(aff1.coeffs)
     n2 = length(aff2.coeffs)
     vars1 = aff1.vars
     vars2 = aff2.vars
     coeffs1 = aff1.coeffs
     coeffs2 = aff2.coeffs
     con1 = aff1.constant
     con2 = aff2.constant

     for i in 1:n1
     	 for j in 1:n2
	     index = -1
	     for k in 1:length(quad.qcoeffs)
	     	 if (quad.qvars1[k] == vars1[i] && quad.qvars2[k] == vars2[j]) || (quad.qvars1[k] == vars2[i] && quad.qvars2[k] == vars1[j])
		     index = k
		     break
		 end
	     end
	     if index != -1
	     	 quad.qcoeffs[k] += coeffs1[i]*coeffs2[j]*c
		 break	     
	     end
	     if coeffs1[i]*coeffs2[j]*c != 0
	         push!(quad.qvars1, copy(vars1[i], vars1[i].m))
	         push!(quad.qvars2, copy(vars2[j], vars2[j].m))
	         push!(quad.qcoeffs, coeffs1[i]*coeffs2[j]*c)
	     end	  
	 end
     end
     temp = addAff(scaleAff(copy(aff1), c*con2), scaleAff(copy(aff2), c*con1))
     
     for i in 1:length(temp.coeffs)
     	 if temp.coeffs[i] != 0.0
	     push!(quad.aff.coeffs, temp.coeffs[i])
	     push!(quad.aff.vars, temp.vars[i])
	 end
     end
     quad.aff.constant = c*con1*con2
     return quad
end

function addQuad(a::QuadExpr, b::QuadExpr)    	 
    a.qcoeffs = [a.qcoeffs; b.qcoeffs]
    a.qvars1 = [a.qvars1; b.qvars1]
    a.qvars2 = [a.qvars2; b.qvars2]
    addAff(a.aff, b.aff)
end

function minusQuad(a::QuadExpr, b::QuadExpr)
    a.qcoeffs = [a.qcoeffs; -b.qcoeffs]
    a.qvars1 = [a.qvars1; b.qvars1]
    a.qvars2 = [a.qvars2; b.qvars2]
    minusAff(a.aff, b.aff)
end



function parselinear(expr)
    @assert islinear(expr) == true
    #aff = AffExpr(Variable[],Float64[],0.0)
    constant = 0
    varId = -1
    coeff = 0

    if isconstant(expr)
       constant += eval(expr)
       return varId, coeff, constant
    end
    if typeof(expr) == Expr
       coeff = 1
       if expr.head == :ref
	   varId = expr.args[2]
       end
       if  expr.args[1] == :(*)
	    coeff = 1
            for i = 2:length(expr.args)
                if  isconstant(expr.args[i])
		    coeff = coeff*eval(expr.args[i])
                    continue
                elseif  typeof(expr.args[i]) == Expr
                    if  expr.args[i].head == :ref
		    	varId = expr.args[i].args[2]
                    end
                end
            end
        end
    end
    return varId, coeff, constant
end


function parseAff(m, expr)
    #println("inside parseAff   ",expr)
    @assert isaff(expr) == true	 
    aff = AffExpr(Variable[],Float64[], 0.0)    
    if islinear(expr)
        varid, coeff, constant = parselinear(expr)
	#=
	aff.constant = constant
	if varid != -1
	    aff.vars = [Variable(m, varid)]
	    aff.coeffs = [coeff]
	end
	=#
	#println(varid, coeff, constant)
	if varid == -1
	    aff = AffExpr(Variable[],Float64[],constant)
	else
	    aff = AffExpr([Variable(m, varid)],[coeff],constant)
	end    
	return aff
    end
    if typeof(expr) == Expr
        if  (expr.args[1] == :(+))
            for i = 2:length(expr.args)
	    	addAff(aff, parseAff(m, expr.args[i]))	       	
            end
	elseif expr.args[1] == :(-)
	     if length(expr.args) == 3
	     	addAff(aff, parseAff(m, expr.args[2]))  
	     	minusAff(aff, parseAff(m, expr.args[3]))    
	     elseif length(expr.args) == 2
	     	minusAff(aff, parseAff(m, expr.args[2]))    
	     end	
        end
    end
    return aff
end

function parseQuad(m, expr)
    #println("inside parseQuad   ",expr)
    @assert isQuad(expr) == true
    quad = QuadExpr()

    if isaff(expr)
       aff = parseAff(m, expr)
       quad.aff = copy(aff)
       return quad
    end
    if typeof(expr) == Expr
       if expr.args[1] == :(*)       	    
       	    #=
	    c = 1
            for i = 2:length(expr.args)
                if isconstant(expr.args[i])
                    c *= eval(expr.args[i])
		elseif islinear(expr.args[i])
		    varid, coeff, ~ = parselinear(expr.args[i])
		    c *= coeff
		    if length(quad.qvars1) == 0
		       push!(quad.qvars1, Variable(m, varid))
		    else
		       push!(quad.qvars2, Variable(m, varid))	
		    end   
                    continue
                else
                    return false
                end
            end
	    =#
	    index = []
            c = 1
            for i = 2:length(expr.args)
                if isconstant(expr.args[i])
                    c *= eval(expr.args[i])
                elseif isaff(expr.args[i])
		    push!(index, i)      
                else
		    error("not quadratic")
                    return false
                end
            end
	    aff1 = parseAff(m, expr.args[index[1]])
	    aff2 = parseAff(m, expr.args[index[2]])	    
	    quad = multAff(m, aff1, aff2, c)
       end
       if  expr.args[1] == :(+) 
            for i = 2:length(expr.args)
	    	addQuad(quad, parseQuad(m, expr.args[i]))
            end
       elseif expr.args[1] == :(-)
             if length(expr.args) == 3
                addQuad(quad, parseQuad(m, expr.args[2]))
                minusQuad(quad, parseQuad(m, expr.args[3]))
             elseif length(expr.args) == 2
                minusQuad(quad, parseQuad(m, expr.args[2]))
             end
    
        end
    end
    #println("quad:    ",quad)
    return quad
end


function _modifycon!(m, expr::Expr)
    #println("inside modifycon!    ", expr)	 
    @assert length(expr.args) == 3
    mainex = expr.args[2]
    sense = expr.args[1]
    if isaff(mainex)
        exprtype = :(aff)
    elseif isQuad(mainex)
	exprtype = :(quad)  
    elseif isexponential(mainex) || islog(mainex) || ispower(mainex) || ismonomial(mainex)	
    	exprtype = :(notaff)   
	_splicevars!(expr, m)
	addNLconstraint2(m, expr)
	return 
    else	 
        exprtype, mainex = _addcon!(m, mainex, :(quad))
    end

    if exprtype == :(aff)
        aff = parseAff(m, mainex)
        if sense == :(<=)
            @constraint(m, sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant <= expr.args[3])
        elseif sense == :(>=)
            @constraint(m, sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant >= expr.args[3])
        else
            @constraint(m, sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant == expr.args[3])
        end
    elseif exprtype == :(quad)
        quad = parseQuad(m, mainex)
        qvars1 = quad.qvars1
        qvars2 = quad.qvars2
        qcoeffs = quad.qcoeffs
        if sense == :(<=)
            @constraint(m, quad.aff + sum(quad.qcoeffs[i]*quad.qvars1[i]*quad.qvars2[i] for i in 1:length(quad.qcoeffs)) <= expr.args[3])
	    elseif sense == :(>=)
            @constraint(m, quad.aff + sum(quad.qcoeffs[i]*quad.qvars1[i]*quad.qvars2[i] for i in 1:length(quad.qcoeffs)) >= expr.args[3])
        else
            @constraint(m, quad.aff + sum(quad.qcoeffs[i]*quad.qvars1[i]*quad.qvars2[i] for i in 1:length(quad.qcoeffs)) == expr.args[3])
        end
    end
end


function _addcon!(m, expr, output = :(aff))
    #println("inside _addcon!  ", expr, "      ",isaff(expr),"     ",  isexponential(expr))
    if islinear(expr)
        return :(aff), expr
    end

    if output == :(quad) && isQuad(expr)
	    return :(quad), expr		
    end

    if (expr.args[1] == :(+) || expr.args[1] == :(-) )
        realoutput = :(aff)
        for i = 2:length(expr.args)
            childexpr = copy(expr.args[i])
	    if output == :(quad)
	        if !isQuad(childexpr)
                    ~, expr.args[i]  = _addcon!(m, expr.args[i], :(quad))
		    if !isaff(childexpr)
		        realoutput = :(quad)
		    end
            	end
	    elseif output == :(aff)
	        if !isaff(childexpr)
	            ~, expr.args[i]  = _addcon!(m, expr.args[i])
                end
            else
		error("wrong output format")		
	    end	
        end
        return realoutput, expr
    elseif expr.args[1] == :(*)
    	#process exp(x[1])*log(x[2])   	
        for i = 2:length(expr.args)
            childexpr = copy(expr.args[i])
            if !isaff(childexpr)
		#println("hello   ", expr.args[i])
                ~, expr.args[i]  = _addcon!(m, expr.args[i])
		#println("hello2   ", expr.args[i])
            end
        end
	# process x[1]*x[2]*x[3]
	# println("hi  ", expr)
        n = 0
        index = [] # index of nonconstant
        for i = 2:length(expr.args)
            if !isconstant(expr.args[i])
	        n = n + 1
		push!(index, i)
		if n >= 3
		    break
		end		
            end
        end
        if output == :(quad) && n ==2
            return :(quad), expr
        end

	if n >= 2
	    ex1 = expr.args[index[1]]
	    ex2 = expr.args[index[2]]

	    if islinear(ex1) && islinear(ex2)
	        @assert islinear(ex1) == true
	        @assert islinear(ex2) == true
	        varid1, coeff1, ~ = parselinear(ex1)
	        varid2, coeff2, ~ = parselinear(ex2)

	        coeff = coeff1*coeff2
	        var1 = Variable(m, varid1)
	        var2 = Variable(m, varid2)
		bilinear = @variable(m, start = m.colVal[varid1]*m.colVal[varid2])
                xl=getlowerbound(var1)
                xu=getupperbound(var1)
                yl=getlowerbound(var2)
                yu=getupperbound(var2)
                m.colLower[end] = min(xl*yl, xl*yu, xu*yl, xu*yu)
                m.colUpper[end] = max(xl*yl, xl*yu, xu*yl, xu*yu)

                numCols = m.numCols
                newvar = :(x[$numCols])
                @constraint(m, bilinear == var1*var2)	
                expr.args[index[1]] = coeff
                expr.args[index[2]] = newvar

	        #@constraint(m, bilinear == coeff*var1*var2) 
	        #expr.args[index[1]] = 1.0
	        #expr.args[index[2]] = newvar
	    else
                bilinear = @variable(m)
                numCols = m.numCols
                newvar = :(x[$numCols])
                newexpr = Expr(:call, :*, copy(ex1), copy(ex2))
		newexpr = Expr(:call, :-, newexpr, newvar)
                newexpr = Expr(:call, :(==), newexpr, 0)
                _modifycon!(m, newexpr)
                expr.args[index[1]] = 1.0
                expr.args[index[2]] = newvar
	    end
	    #println(m.quadconstr[end])
	end	
	if n >= 3	
            ~, expr  = _addcon!(m, expr)
	end    
	#println("end of *", expr)
        return :(aff), expr

    elseif expr.args[1] == :(/)
              new = @variable(m)
              numCols = m.numCols
              newvar = :(x[$numCols])
              newexpr = Expr(:call, :*, copy(expr.args[3]), newvar)
	      newexpr = Expr(:call, :-, newexpr, copy(expr.args[2]))
              newexpr = Expr(:call, :(==), newexpr, 0)
              _modifycon!(m, newexpr)
              expr = newvar
              return :(aff), expr

    elseif expr.args[1] == :(exp) 
	      if !islinear(expr.args[2])	      
		  new = @variable(m)
		  numCols = m.numCols
		  newvar = :(x[$numCols])
		  
 		  newexpr = Expr(:call, :-, copy(expr.args[2]), newvar)
              	  newexpr = Expr(:call, :(==), newexpr, 0)
		  _modifycon!(m, newexpr)
		  expr.args[2] = newvar
	      end	      		
              new = @variable(m)
              numCols = m.numCols
              newvar = :(x[$numCols])
              newexpr = Expr(:call, :-, copy(expr), newvar)
              newexpr = Expr(:call, :(==), newexpr, 0)
              _modifycon!(m, newexpr)
	      expr = newvar
	      return :(exp), expr	      

    elseif expr.args[1] == :(log)
	      #println("inside log   ", expr)
              if !islinear(expr.args[2])
                  new = @variable(m)
                  numCols = m.numCols
                  newvar = :(x[$numCols])
                  newexpr = Expr(:call, :-, copy(expr.args[2]), newvar)
                  newexpr = Expr(:call, :(==), newexpr, 0)
                  _modifycon!(m, newexpr)
                  expr.args[2] = newvar
              end

              new = @variable(m)
              numCols = m.numCols
              newvar = :(x[$numCols])
              newexpr = Expr(:call, :-, copy(expr), newvar)
              newexpr = Expr(:call, :(==), newexpr, 0)
              _modifycon!(m, newexpr)
              expr = newvar
	      #println("inside log   ", expr)
              return :(aff), expr

    elseif expr.args[1] == :(^) 
    	      if isconstant(expr.args[2])	  #a^(bx[1])    	 
              	  if !islinear(expr.args[3])	
                      new = @variable(m)
                      numCols = m.numCols
                      newvar = :(x[$numCols])	     
                      newexpr = Expr(:call, :-, copy(expr.args[3]), newvar)
                      newexpr = Expr(:call, :(==), newexpr, 0)
                      _modifycon!(m, newexpr)
                      expr.args[3] = newvar
		  end 
	      elseif isconstant(expr.args[3])      #(bx[1])^a
	          d = expr.args[3]
                  if output == :(quad) && isaff(expr.args[2]) && d == 2
                            expr.args[1] = :(*)
		            expr.args[3] = copy(expr.args[2])			
                            return :(quad), expr
                  end
                  if !islinear(expr.args[2])
                      new = @variable(m)
                      numCols = m.numCols
                      newvar = :(x[$numCols])
                      newexpr = Expr(:call, :-, copy(expr.args[2]), newvar)
                      newexpr = Expr(:call, :(==), newexpr, 0)
                      _modifycon!(m, newexpr)
                      expr.args[2] = newvar
                  end	     
              end
	      
              if output == :(quad) &&   isaff(expr.args[2]) && d == 2
	      	  expr.args[1] = :(*)
		  expr.args[3] = copy(expr.args[2])	      	                    
		  return :(quad), expr				
              end

              new = @variable(m)
              numCols = m.numCols
              newvar = :(x[$numCols])
              newexpr = Expr(:call, :-, copy(expr), newvar)
              newexpr = Expr(:call, :(==), newexpr, 0)
              _modifycon!(m, newexpr)
              expr = newvar
              return :(aff), expr
    else
        return :(error)
    end
end

#:(x[1] - exp(b*x[2])) 
function isexponential(expr::Expr)
    if typeof(expr) == Expr

        if length(expr.args) == 3 && (expr.args[1] == :(-)) # || expr.args[1] == :(+))
            if islinear(expr.args[2]) && isconstant(expr.args[2]) == false
	        lloc = 2
		linear = expr.args[2]
	        left = expr.args[3]
            elseif islinear(expr.args[3]) &&  (isconstant(expr.args[3]) == false)
	    	lloc =3
		linear = expr.args[3]   
		left = expr.args[2]
	    else
		return false	
	    end
	    lvarId, a, constant = parselinear(linear)
	    if a != 1.0 || constant!=0
	        return false
	    end

    	    if typeof(left) == Expr
	        c = 1
		n = 0
       	        if left.args[1] == :exp && islinear(left.args[2]) && isconstant(left.args[2])==false   #left.args[2].head == :ref
		    n = 1
                    temp = left.args[2]
                end
		if left.args[1] == :*
            	    for i = 2:length(left.args)
                    	if  isconstant(left.args[i])
                    	    c = c*eval(left.args[i])
                    	    continue
                	elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :exp && islinear(left.args[i].args[2]) && isconstant(left.args[i].args[2]) == false
			    temp = left.args[i].args[2]	
			    n += 1
			else
			    return false      
                        end
                    end
		end    
		if n == 1 && c == 1.0
		        lvarId, b, constant = parselinear(temp)
			if b != 0 && constant == 0 
		            return true
			end    
                end
            end
        end
    end
    return false
end


#:(x[1] - log(b*x[2]))
function islog(expr::Expr)
    if typeof(expr) == Expr
        if length(expr.args) == 3 && (expr.args[1] == :(-)) # || expr.args[1] == :(+))
            if islinear(expr.args[2]) && isconstant(expr.args[2]) == false
                lloc = 2
                linear = expr.args[2]
                left = expr.args[3]
            elseif islinear(expr.args[3]) &&  (isconstant(expr.args[3]) == false)
                lloc =3
                linear = expr.args[3]
                left = expr.args[2]
            else
                return false
            end
            lvarId, a, constant = parselinear(linear)
            if a != 1.0 || constant!=0
                return false
            end
            if typeof(left) == Expr
                c = 1
                n = 0
                if left.args[1] == :log && islinear(left.args[2]) && isconstant(left.args[2])==false   #left.args[2].head == :ref
                    n = 1
                    temp = left.args[2]
                end
                if left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :log && islinear(left.args[i].args[2]) && isconstant(left.args[i].args[2]) == false
                            temp = left.args[i].args[2]
                            n += 1
                        else
                            return false
                        end
                    end
                end
                if n == 1 && c == 1.0
                        lvarId, b, constant = parselinear(temp)
                        if b != 0 && constant == 0
                            return true
                        end
                end
            end
        end
    end
    return false
end


#:(x[1] - (b*x[2])^d)
function ispower(expr::Expr)
    if typeof(expr) == Expr
        if length(expr.args) == 3 && (expr.args[1] == :(-)) # || expr.args[1] == :(+))
            if islinear(expr.args[2]) && isconstant(expr.args[2]) == false
                lloc = 2
                linear = expr.args[2]
                left = expr.args[3]
            elseif islinear(expr.args[3]) &&  (isconstant(expr.args[3]) == false)
                lloc =3
                linear = expr.args[3]
                left = expr.args[2]
            else
                return false
            end
            lvarId, a, constant = parselinear(linear)
            if a != 1.0 || constant!=0
                return false
            end

            if typeof(left) == Expr
                c = 1
                n = 0
                if left.args[1] == :^ && islinear(left.args[2]) && isconstant(left.args[3])   #left.args[2].head == :ref
                    n = 1
		    left.args[3] = eval(left.args[3])
                    temp = left.args[2]
                end
                if left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :^ && islinear(left.args[i].args[2]) && isconstant(left.args[i].args[3])
			    left.args[i].args[3] = eval(left.args[i].args[3])	
                            temp = left.args[i].args[2]
                            n += 1
                        else
                            return false
                        end
                    end
                end
                if n == 1 && c == 1.0
                        lvarId, b, constant = parselinear(temp)
                        if b != 0 && constant == 0
                            return true
                        end
                end
            end
        end
    end
    return false
end


#:(x[1] - (d)^(b*x[2]))
function ismonomial(expr::Expr)
    if typeof(expr) == Expr
        if length(expr.args) == 3 && (expr.args[1] == :(-)) # || expr.args[1] == :(+))
            if islinear(expr.args[2]) && isconstant(expr.args[2]) == false
                lloc = 2
                linear = expr.args[2]
                left = expr.args[3]
            elseif islinear(expr.args[3]) &&  (isconstant(expr.args[3]) == false)
                lloc = 3
                linear = expr.args[3]
                left = expr.args[2]
            else
                return false
            end
            lvarId, a, constant = parselinear(linear)
            if a != 1.0 || constant!=0
                return false
            end

            if typeof(left) == Expr
                c = 1
                n = 0
                if left.args[1] == :^ && islinear(left.args[3]) && isconstant(left.args[2])   #left.args[2].head == :ref
                    n = 1
		    left.args[2] = eval(left.args[2])  
                    temp = left.args[3]
                end
                if left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :^ && islinear(left.args[i].args[3]) && isconstant(left.args[i].args[2])
			    left.args[i].args[3] = eval(left.args[i].args[3])	
                            temp = left.args[i].args[3]
                            n += 1
                        else
                            return false
                        end
                    end
                end
                if n == 1 && c == 1.0
                        lvarId, b, constant = parselinear(temp)
                        if b != 0 && constant == 0
                            return true
                        end
                end
            end
        end
    end
    return false
end


# :(x[1] - exp(x[2]) <= 0.0)
# NLconstraint(m, x[1] <= exp(x[2]))
function isexponentialcon(expr::Expr)
    @assert length(expr.args) == 3
    mainex = expr.args[2]
    return isexponential(mainex)
end

function islogcon(expr::Expr)
    @assert length(expr.args) == 3
    mainex = expr.args[2]
    return islog(mainex)
end

function ispowercon(expr::Expr)
    @assert length(expr.args) == 3
    mainex = expr.args[2]
    return ispower(mainex)
end

function ismonomialcon(expr::Expr)
    @assert length(expr.args) == 3
    mainex = expr.args[2]
    return ismonomial(mainex)
end


#(x[1] op exp(b*x[2]))
# assume already check by isexponential
function parseexponentialcon(expr::Expr)
    @assert length(expr.args) == 3
    @assert expr.args[3] == 0.0
    mainex = expr.args[2]
    @assert mainex.args[1] == :(-)

    if islinear(mainex.args[2])
       lloc = 2
       linear = mainex.args[2]
       left = mainex.args[3]
    else
       lloc = 3	
       linear = mainex.args[3]
       left = mainex.args[2]
    end          
    lvarId, a, c = parselinear(linear)
    @assert a == 1.0    
    @assert c == 0.0

    c = 1
    if left.args[1] == :exp && islinear(left.args[2])
         temp = left.args[2]
    elseif left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :exp && islinear(left.args[i].args[2])
                            temp = left.args[i].args[2]
			end
                    end
    end
    @assert c == 1.0
    nlvarId, b, ~ = parselinear(temp)	
    if expr.args[1] == :(==)
           op = :(==)
    elseif expr.args[1] == :(<=)
           if lloc == 2
              op = :(<=)
           else
              op = :(>=)
           end
    elseif expr.args[1] == :(>=)
           if lloc == 2
               op = :(>=)
           else
               op = :(<=)
           end
    end
    return ExpVariable(lvarId, nlvarId, op, b, [])
end


#(x[1] op log(b*x[2]))
# assume already check by isexponential
function parselogcon(expr::Expr)
    @assert length(expr.args) == 3
    @assert expr.args[3] == 0.0
    mainex = expr.args[2]
    @assert mainex.args[1] == :(-)

    if islinear(mainex.args[2])
       lloc = 2
       linear = mainex.args[2]
       left = mainex.args[3]
    else
       lloc = 3
       linear = mainex.args[3]
       left = mainex.args[2]
    end
    lvarId, a, c = parselinear(linear)
    @assert a == 1.0
    @assert c == 0.0

    c = 1
    if left.args[1] == :log && islinear(left.args[2])
         temp = left.args[2]
    elseif left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :log && islinear(left.args[i].args[2])
                            temp = left.args[i].args[2]
                        end
                    end
    end
    @assert c == 1.0
    nlvarId, b, ~ = parselinear(temp)
    if expr.args[1] == :(==)
           op = :(==)
    elseif expr.args[1] == :(<=)
           if lloc == 2
              op = :(<=)
           else
              op = :(>=)
           end
    elseif expr.args[1] == :(>=)
           if lloc == 2
               op = :(>=)
           else
               op = :(<=)
           end
    end
    return LogVariable(lvarId, nlvarId, op, b, [])
end


#:(x[1] op (b*x[2])^d)
function parsepowercon(expr::Expr)
    #println("inside parsepower ", expr)
    @assert length(expr.args) == 3
    @assert expr.args[3] == 0.0
    mainex = expr.args[2]
    @assert mainex.args[1] == :(-)

    if islinear(mainex.args[2])
       lloc = 2
       linear = mainex.args[2]
       left = mainex.args[3]
    else
       lloc = 3
       linear = mainex.args[3]
       left = mainex.args[2]
    end
    lvarId, a, c = parselinear(linear)
    @assert a == 1.0
    @assert c == 0.0

    c = 1
    if left.args[1] == :^ && islinear(left.args[2]) && isconstant(left.args[3]) 
         d = eval(left.args[3])
         temp = left.args[2]
    elseif left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :^ && islinear(left.args[i].args[2]) && isconstant(left.args[i].args[3])
			    d = eval(left.args[i].args[3])	
                            temp = left.args[i].args[2]
			end
                    end
    end
    @assert c == 1.0
    nlvarId, b, ~ = parselinear(temp)

    if expr.args[1] == :(==)
           op = :(==)
    elseif expr.args[1] == :(<=)
           if lloc == 2
              op = :(<=)
           else
              op = :(>=)
           end
    elseif expr.args[1] == :(>=)
           if lloc == 2
               op = :(>=)
           else
               op = :(<=)
           end
    end
    #println(lvarId, "   ", nlvarId, "    ",op, "   ",b, "     ",d)
    return PowerVariable(lvarId, nlvarId, op, b, d, [])
end



#:(x[1] op d^(b*x[2])
function parsemonomialcon(expr::Expr)
    @assert length(expr.args) == 3
    @assert expr.args[3] == 0.0
    mainex = expr.args[2]
    @assert mainex.args[1] == :(-)

    if islinear(mainex.args[2])
       lloc = 2
       linear = mainex.args[2]
       left = mainex.args[3]
    else
       lloc = 3
       linear = mainex.args[3]
       left = mainex.args[2]
    end
    lvarId, a, c = parselinear(linear)
    @assert a == 1.0
    @assert c == 0.0

    c = 1
    if left.args[1] == :^ && islinear(left.args[3]) && isconstant(left.args[2])
         d = eval(left.args[2])
         temp = left.args[3]
    elseif left.args[1] == :*
                    for i = 2:length(left.args)
                        if  isconstant(left.args[i])
                            c = c*eval(left.args[i])
                            continue
                        elseif  typeof(left.args[i]) == Expr && left.args[i].args[1] == :^ && islinear(left.args[i].args[3]) && isconstant(left.args[i].args[2])
                            d = eval(left.args[i].args[2])
                            temp = left.args[i].args[3]
                        end
                    end
    end
    @assert c == 1.0
    nlvarId, b, ~ = parselinear(temp)

    if expr.args[1] == :(==)
           op = :(==)
    elseif expr.args[1] == :(<=)
           if lloc == 2
              op = :(<=)
           else
              op = :(>=)
           end
    elseif expr.args[1] == :(>=)
           if lloc == 2
               op = :(>=)
           else
               op = :(<=)
           end
    end


    if d <= 0
         error("warning:  a^x, a cannot be negative, please reformulate ")
    end
    b = b*log(d)
    return ExpVariable(lvarId, nlvarId, op, b, [])     
    #return MonomialVariable(lvarId, nlvarId, op, b, d, [])
end



function isanumber(x)
    if typeof(x) == Int || typeof(x) == Float64
       return true
    end
    return false
end

# true if "1+1"
#x not comparison
function isconstant(x) 
    if isanumber(x)
       return true
    end
    if typeof(x) == Expr
       if x.head == :ref
       	   return false	   
       end
       for i = 1:length(x.args)
       	   if typeof(x.args[i]) == Expr
	      if x.args[i].head == :ref  
	      	 return false
	      else 
	      	 if ! isconstant(x.args[i])
		    return false
		 end	 
              end
           end
       end
       return true
    end
    return false
end



# true if "2*x[1]*x[3]+5*x[1]"
function isQuad(expr)
    #println("inside isquad ", expr)
    if isaff(expr)
       return true
    end
    if typeof(expr) == Expr
       if expr.args[1] == :(*) 
            n = 0
	    #println("inside *")
            for i = 2:length(expr.args)
	    	#println("args[i]:  ",expr.args[i])	    	
	    	if isconstant(expr.args[i])
		    continue
                elseif isaff(expr.args[i])
                    n += 1
                    continue
                else
                    return false
                end
		#=    
                elseif islinear(expr.args[i]) 
		    ~, ~, constant = parselinear(expr.args[i])
		    if constant != 0
		       	return false
		    end
		    n += 1
                    continue
                else
                    return false
                end
		=#
            end
	    #println("n",n)
	    if n <= 2
	       return true
	    end
       end
       if  expr.args[1] == :(+) || expr.args[1] == :(-)
            for i = 2:length(expr.args)
                if  isaff(expr.args[i]) || isQuad(expr.args[i])
		    continue	
		else
		    return false    
                end
            end
	    #println("is quad")
	    return true
        end
    end
    return false
end


# true if "5*x[1]"
function islinear(expr)
    if isconstant(expr)
       return true
    end
    if typeof(expr) == Expr
       if expr.head == :ref
           return true
       end
	if  expr.args[1] == :(*) 
    	    n = 0
    	    for i = 2:length(expr.args)
	    	if  isconstant(expr.args[i])
		    continue
	    	elseif  typeof(expr.args[i]) == Expr
	       	    if  expr.args[i].head == :ref
	       	    	n +=1
		    else
			return false	
                    end
                end
	    end 
            if  n <= 1
                return true
            end 
	end
    end
    return false
end

# true if "5*x[1]+6*x[2]+7"
function isaff(expr)
    #println("inside isaff  ",expr)	 
    if islinear(expr)
       return true
    end
    if typeof(expr) == Expr
        if  (expr.args[1] == :(+) || expr.args[1] == :(-))
            for i = 2:length(expr.args)
                if  islinear(expr.args[i])
                    continue
                elseif isaff(expr.args[i])
	            continue
		else
		    return false 
                end
            end
 	    return true
        end

        if  (expr.args[1] == :(*))
            n = 0	
            for i = 2:length(expr.args)
                if  isconstant(expr.args[i])
                    continue
                elseif isaff(expr.args[i])
		    n += 1
                    continue
		    
                else
                    return false
                end
            end
            if  n <= 1
                return true
            end
        end
    end
    return false
end	      


function hasBinaryVar(m::JuMP.Model)
    for i = 1:m.numCols
    	if m.colCat[i] == :Bin
    	    return true
	end
    end	
    return false
end


function numBinaryVar(m::JuMP.Model)
    n = 0
    for i = 1:m.numCols
        if m.colCat[i] == :Bin
            n = n + 1
        end
    end
    return n
end


function fixBinaryVar(m::JuMP.Model)
    for i = 1:m.numCols
        if m.colCat[i] == :Bin
            val = 0.0
            if m.colVal[i] >= 0.5
                val = 1.0
            end
            m.colCat[i] = :Fixed
            m.colLower[i] = val
            m.colUpper[i] = val
            m.colVal[i] = val
        end
    end
    
    for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
    	for i = 1:scenario.numCols
            if scenario.colCat[i] == :Bin
                val = 0.0
            	if scenario.colVal[i] >= 0.5
                    val = 1.0
                end
            	scenario.colCat[i] = :Fixed
            	scenario.colLower[i] = val
            	scenario.colUpper[i] = val
            	scenario.colVal[i] = val
            end
        end
    end
end
#=
function fixBinaryVar(m::JuMP.Model, sol)
    for i = 1:m.numCols
    if m.colCat[i] == :Bin
            val = 0.0
            if sol[i] >= 0.5
	       val = 1.0
            end
            m.colCat[i] = :Fixed
            m.colLower[i] = val
            m.colUpper[i] = val
            m.colVal[i] = val
	    end
    end
end
=#
function fixVar(v::Variable, val::Number)
    v.m.colCat[v.col] = :Fixed
    v.m.colLower[v.col] = val
    v.m.colUpper[v.col] = val
    v.m.colVal[v.col] = val
end

#=
function factorable(P)
   return P
end
=#
include("boundT.jl")
include("relax.jl")
include("updaterelax.jl")
include("preprocessex.jl")
include("preprocessSto.jl")
include("preprocess.jl")