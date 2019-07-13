const n_rel = 100
const lambda = 0.25
const alpha = 0.1
const mu_score = 0.15
const lambda_consecutive = 4
const sigma_violation = 1e-8
const con_tol = 1e-6
const LP_improve_tol = 1e-2
const mingap = 1e-2

const machine_error = 1e-10
const small_bound_improve = 1e-4
const large_bound_improve = 1e-1
const probing_improve = 1e-4
const nR = 0
const probing = true
const max_score_min = 1e-3

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
end
Node() = Node(Float64[], Float64[], -1, -1, -1e10, 0, nothing, nothing, 0)

type PreprocessResult
     branchVarsId
     EqVconstr
     multiVariable_list
     multiVariable_convex
     multiVariable_aBB
     qbVarsId
end
PreprocessResult() = PreprocessResult(Int[], LinearConstraint[], [], [], [], [])


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
            if firstVarsId[i] > 0
                if scenario.colLower[firstVarsId[i]] >  P.colLower[i]
                    P.colLower[i] = scenario.colLower[firstVarsId[i]]
                end
                if scenario.colUpper[firstVarsId[i]] <  P.colUpper[i]
                    P.colUpper[i] = scenario.colUpper[firstVarsId[i]]
                end
            end
        end
    end
    updateFirstBounds!(P, P.colLower, P.colUpper)
end

function updateFirstBounds!(P, xl, xu)
        P.colLower = copy(xl)
        P.colUpper = copy(xu)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))	
	    for j = 1:length(scenario.ext[:firstVarsId])			
	    	if scenario.ext[:firstVarsId][j] != -1
            	    scenario.colLower[scenario.ext[:firstVarsId][j]] = copy(xl[j])
           	    scenario.colUpper[scenario.ext[:firstVarsId][j]] = copy(xu[j])
                end
	    end
        end
end
function updateFirstBounds!(P, xl, xu, varId)
        P.colLower[varId] = copy(xl)
        P.colUpper[varId] = copy(xu)
	for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
	    if scenario.ext[:firstVarsId][varId] != -1
                scenario.colLower[scenario.ext[:firstVarsId][varId]] = copy(xl)
                scenario.colUpper[scenario.ext[:firstVarsId][varId]] = copy(xu)
	    end   
        end
end


function updateExtensiveFirstBounds!(Pex, P, xl, xu)
        n = P.numCols
        Pex.colLower[1:n] = copy(xl)
        Pex.colUpper[1:n] = copy(xu)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            for j = 1:length(scenario.ext[:firstVarsId])
                if scenario.ext[:firstVarsId][j] != -1
                    Pex.colLower[v_map[scenario.ext[:firstVarsId][j]]] = copy(xl[j])
                    Pex.colUpper[v_map[scenario.ext[:firstVarsId][j]]] = copy(xu[j])
                end
            end
            n = n + scenario.numCols
        end
end


function updateExtensiveFirstBounds!(Pex, P, xl, xu, varId)
        n = P.numCols
        Pex.colLower[varId] = copy(xl)
        Pex.colUpper[varId] = copy(xu)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
	    if scenario.ext[:firstVarsId][varId] != -1
                Pex.colLower[v_map[scenario.ext[:firstVarsId][varId]]] = copy(xl)
            	Pex.colUpper[v_map[scenario.ext[:firstVarsId][varId]]] = copy(xu)
	    end	
        end
end

function updateExtensiveFirstBounds!(PexLower, PexUpper, P, Pex, xl, xu, varId)
        n = P.numCols
        PexLower[varId] = copy(xl)
        PexUpper[varId] = copy(xu)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
	    if scenario.ext[:firstVarsId][varId] != -1
                PexLower[v_map[scenario.ext[:firstVarsId][varId]]] = copy(xl)
                PexUpper[v_map[scenario.ext[:firstVarsId][varId]]] = copy(xu)
            end
        end
end

function updateExtensiveBoundsFromSto!(P, Pex)
        nfirst = P.numCols
        Pex.colLower[1:nfirst] = copy(P.colLower)
       	Pex.colUpper[1:nfirst] = copy(P.colUpper)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            Pex.colLower[v_map] = copy(scenario.colLower)
            Pex.colUpper[v_map] = copy(scenario.colUpper)
        end
end

function updateStoBoundsFromExtensive!(Pex, P)
	updateStoFirstBounds!(P) 
        nfirst = P.numCols
        P.colLower = copy(Pex.colLower[1:nfirst])
        P.colUpper = copy(Pex.colUpper[1:nfirst])
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(P))
            v_map = Pex.ext[:v_map][idx]
            scenario.colLower = copy(Pex.colLower[v_map])
            scenario.colUpper = copy(Pex.colUpper[v_map])
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
	mchildrens = PlasmoOld.getchildren(m)
        for (idx,scenario) in enumerate(PlasmoOld.getchildren(dest))
            scenario.colLower = copy(mchildrens[idx].colLower)
            scenario.colUpper = copy(mchildrens[idx].colUpper)
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

include("boundT.jl")
include("relax.jl")
include("updaterelax.jl")
include("preprocess.jl")
include("preprocessex.jl")
include("preprocessSto.jl")







#=
function decompose(P)
    #println("relax1")
    #m = Model(solver=IpoptSolver(print_level = 0))
    #m = Model(solver = GLPKSolverLP())
    #m.solver = P.solver  # The two models are linked by this

    # Variables
    m.numCols = P.numCols
    m.colNames = P.colNames[:]
    m.colNamesIJulia = P.colNamesIJulia[:]
    m.colLower = P.colLower[:]
    m.colUpper = P.colUpper[:]
    #m.colLower = copy(xl)
    #m.colUpper = copy(xu)
    m.colCat = P.colCat[:]
    m.colVal = P.colVal[:]
    m.varDict = Dict{Symbol,Any}()
    for (symb,v) in P.varDict
        m.varDict[symb] = copy(v, m)
    end
    # Constraints
    m.linconstr  = map(c->copy(c, m), P.linconstr)
    m.quadconstr = map(c->copy(c, m), P.quadconstr)

    # Objective
    m.obj = copy(P.obj, m)
    m.objSense = P.objSense

      
    if P.nlpdata != nothing
       d = JuMP.NLPEvaluator(P)         #Get the NLP evaluator object.  Initialize the expression graph
       MathProgBase.initialize(d,[:ExprGraph])
       num_cons = MathProgBase.numconstr(node_model)
       for i = (1+length(m.linconstr)+length(m.quadconstr)):num_cons
        expr = MathProgBase.constr_expr(d,i)  #this returns a julia expression
        _modifycon!(m, expr)            #splice the variables from v_map into the expression
       end
    end
    return m
end


function _modifycon!(m, expr::Expr, numcols)	 
    if length(expr.args) == 3	 
       mainex = expr[1]
    elseif length(expr.args) == 5
       mainex = expr[3]		
    end

    _addcon!(m, mainex)


    _splicevars!(m, expr::Expr)

    JuMP.addNLconstraint(m,expr)


    #=
    aff = exprtoaff(mainex)
    if length(expr.args) == 3
       sense = expr.args[2]
       if sense == :(<=)
       	  @constraint(m, sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant <= expr.args[3])
       elseif expr.args[2] == :(>=)
          @constraint(m, sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant >= expr.args[3])
       else
          @constraint(m, sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant == expr.args[3])	 
       end	  
    elseif length(expr.args) == 5
       @constraint(m, expr.args[1] <= sum(aff.coeffs[i]*aff.vars[i] for i = 1:length(aff.coeffs)) + aff.constant <= expr.args[3])   
    end
    =# 
end


function _addcon!(m, expr, numcols)
    if isaff(expr)	 
        return :(aff)      
    elseif isnumber(isexpr.args[2])
    	return :(exp)
    elseif islog(exp)
	return :(log)
    elseif isconstant
    	return :(constant)	
    end 

    if (expr.args[1] == :(+) || expr.args[1] == :(-)) 
    	for i = 2:length(expr.args)
	     childexpr = copy(expr.args[i])
	     if ! isaff(childexpr)   	   
	     	 numcols +=1
		 newvar = :(x[$numcols])
		 expr.args[i] = newvar
		 childexpr = Expr(:call, :-, childexpr, newvar)
		 newexpr = Expr(:comparison, childexpr, :(==), 0)
		 _modifycon!(m, newexpr, numcols) 	     	  
	      end
        end
    end
end


function isanumber(x)
    if typeof(x) == Int || typeof(x) == Float64
       return true
    end
    return false
end

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


function islinear(expr::Expr)
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

function isaff(expr::Expr)
    if islinear(expr)
       return true
    end
    if typeof(expr) == Expr
        if  (expr.args[1] == :(+) || expr.args[1] == :(-))
            for i = 2:length(expr.args)
                if  islinear(expr.args[i])
                    continue
                elseif typeof(expr.args[i]) == Expr
                    if !isaff(expr.args[i]) 
		        return false 
		    end
                end
            end
 	    return true
        end
    end
    return false
end	      

function isexponential(expr::Expr)
    if typeof(expr.args[i]) == Expr	 
       if expr.args[1] == :exp && expr.args[2].head == :ref
       	  return true
       end
    end	 
    return false
end

function islog(expr::Expr)
    if typeof(expr.args[i]) == Expr	
       if expr.args[1] = :log && expr.args[2].head == :ref
          return true
       end
    end
    return false
end

function ispower(expr::Expr)
    if expr.args[1] = :^ && expr.args[2].head != :ref && typeof(expr.args[3]) == Float64

    end
end

function ismonomial(expr::Expr)
    if expr.args[1] = :^ && expr.args[3].head != :ref && typeof(expr.args[2]) == Float64

    end
end

    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splicevars!(expr.args[i],v_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]   #this is the actual index in x[1], x[2], etc...
                new_var = :($(v_map[var_index]))   #get the JuMP variable from v_map using the index
                expr.args[i] = new_var             #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end

   
	 JuMP.addNLconstraint(m,expr) 


end


function _splicevars!(m, expr::Expr)
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splicevars!(expr.args[i],v_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]   #this is the actual index in x[1], x[2], etc...
		new_var = Variable(m, var_index)		
                new_var = :($(new_var))		   #get the JuMP variable from v_map using the index
                expr.args[i] = new_var             #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end
=#



