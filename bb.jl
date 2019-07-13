push!(LOAD_PATH, pwd())
using Gurobi, JuMP, Ipopt, SCIP, PlasmoOld
importall MathProgBase.SolverInterface
using Distributions
EnableNLPResolve()
include("core.jl")


function branch_bound(m)  
    Ipopt_solve(m) 	 	
    P = copyStoModel(m)	  
    scenarios = PlasmoOld.getchildren(P)
    println("No. of first stage variables  ", P.numCols)
    println("No. of second stage variables ", scenarios[1].numCols)
    println("No. of second stage constraints ", MathProgBase.numconstr(scenarios[1]))
    println("No. of scenarios  ", length(scenarios))
    Pcount = extensiveSimplifiedModel(P)
    println("No. of total variables ", Pcount.numCols)
    println("No. of total constraints ", MathProgBase.numconstr(Pcount))
    Pcount = 1

    tic()
    preprocessSto!(P)
    println("preprocessSto:   ",  toq(), " (s)")

    tic()
    scenarios = PlasmoOld.getchildren(P)
    nscen = length(scenarios)
    pr_children = []
    for (idx,scen) in enumerate(scenarios)
        pr = preprocess!(scen)
        push!(pr_children, pr)
    end
    println("preprocess:   ",  toq(), " (s)")
   

    ncols_first = P.numCols
    bVarsId = collect(1:ncols_first)
    vs = VarSelector(bVarsId)

    # P might be in the QP form
    Pprobing = copyStoModel(P)
    P_child = copyStoModel(P)
    PWS = copyStoModel(P)         # P is used for get upper bound from local solver, PWS is used for get lower bound from wait and see solution, PWS is nonlinear form
    PWSfix = copyStoModel(PWS)    # PWSfix is used for get upper bound from wait and see solution with first stage fixed, PWSfixed is nonlinear form
    #PWS = copyNLStoModel(P)	    # P is used for get upper bound from local solver, PWS is used for get lower bound from wait and see solution, PWS is nonlinear form
    #PWSfix = copyNLStoModel(PWS)    # PWSfix is used for get upper bound from wait and see solution with first stage fixed, PWSfixed is nonlinear form
    scenariosWS = PlasmoOld.getchildren(PWS)
    scenariosWSfix = PlasmoOld.getchildren(PWSfix)
    tic()
    Pex = extensiveSimplifiedModel(P)
    Pex_probing = copyModel(Pex)
    println("extensive:   ",  toq(), " (s)")
    #println(Pex)

    tic()
    prex = preprocessex!(Pex)
    println("preprocessex:   ",  toq(), " (s)")

    UB = 1e10
    Pex.solver = IpoptSolver()
    local_status = solve(Pex)
    println(local_status)
    if local_status == :Optimal
        updateStoSolFromExtensive!(Pex, P)
        UB = Pex.objVal
    end
    println("U:   ", UB)


    root = Node(copy(Pex.colLower), copy(Pex.colUpper), -1, 1, -1e10, 0, nothing, nothing, 1e10)
    nodeList =[]
    push!(nodeList, root)

    x = copy(P.colVal)
    WS_status = :Optimal
    relaxed_status = :Optimal
    WSfirst = copy(P.colVal)
    WSSol = StoSol(nscen)

    nfirst = P.numCols
    iter = 0
    println(" iter ", " left ", " lev  ", " bVarId ","      bvlb      ", "       bvub      ", "       LB       ", "       UB      ", "      gap   ")

    tic()
    R=relax(Pex, prex, 1e10) #UB)
    Roriginal = copyModel(R)
    Rold = copyModel(R)
    println("relax:   ",  toq(), " (s)")
    PWS_child = copyStoModel(PWS)
    Pex_child = copyModel(Pex)

    LB = -1e10
    FLB = UB			#lower bound for the fathoned node( with lb <= ub but (ub-lb)<=mingap)
    LB_UB_list = []

    tic()
    while nodeList!=[]    	 
	if iter == 5000
	   break
	end
        iter += 1
	push!(LB_UB_list, [LB, UB, 1, 1])
        if (iter != 1)
            LB_UB_list[iter-1][1] = LB
            LB_UB_list[iter-1][2] = UB
        end

        for (idx,n) in enumerate(nodeList)
            println("remaining ", idx,  "   ", n.LB)
        end
	LB, nodeid = getGlobalLowerBound(nodeList, FLB, UB)
        node = nodeList[nodeid]
	deleteat!(nodeList, nodeid)
        Pex.colLower = copy(node.xl)
        Pex.colUpper = copy(node.xu)
	Rold = copyModel(Roriginal)
	#Rold = Roriginal

	#=
	if node.parRSol != nothing
	    println("  node.parRSol:    ", )
	    Rold.colVal = copy(node.parRSol)
	end
	=#
	updateStoBoundsFromExtensive!(Pex, P)		#updateFirstBounds!(P, node.xl[1:nfirst], node.xu[1:nfirst]) 
	level = node.level
		
        if node.bVarId != -1
           bVarId = node.bVarId
	   @printf "%-6d %-6d %-6d %-6d %-14.4f %-14.4e %-14.4e %-14.4e %-7.4f %s\n" iter length(nodeList) node.level  bVarId P.colLower[bVarId] P.colUpper[bVarId] LB UB (UB-LB)/abs(LB)*100 "%" 	    
        else
	   @printf "%-6d %-6d %-40d %-14.4f %-14.4e %-7.4f %s \n" iter length(nodeList) node.level LB UB (UB-LB)/abs(LB)*100 "%"
        end

        reduction_relax = 1
        node_LB = node.LB
	relaxed_status = :Optimal
	UB_status = :Optimal
	delete_nodes = []	
	feasible = true
        while reduction_relax >= 0.1
            reduction_relax = 0	    	    	    
	    println("before feasibility reduction ", P.colLower, "   ", P.colUpper)
	    tic()
	    feasible = Sto_medium_feasibility_reduction(P, pr_children, Pex, prex, Rold, UB, LB, bVarsId)	
	    #feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, LB)
	    println("after initial feasibility reduction ", P.colLower, "   ", P.colUpper, "  ","   time:   ",  toq(), " (s)")
	    updateExtensiveBoundsFromSto!(P, Pex)

	   
	    if probing && feasible && level <= 1	      
	        tic() 
	        #probing
		n_reduced_probing = 1		 
		#while n_reduced_probing > 0
		i = 1
	    	n_reduced_probing = 0
	    	while i<= nfirst
                    xl = P.colLower[i]
                    xu = P.colUpper[i]
		    println(i,"    ", xl, "    ",xu)
                    updateExtensiveBoundsFromSto!(P, Pex_probing)
                    Rprobing = updaterelax(Rold, Pex_probing, prex, 1e10)
                    Rprobing.colCat[i] = :Fixed
                    Rprobing.colLower[i] = xl
                    Rprobing.colUpper[i] = xl
                    Rprobing.colVal[i] = xl
                    relaxed_status = solve(Rprobing)
                    relaxed_LB = getRobjective(Rprobing, relaxed_status, LB)
                    println("relaxed_status  ", relaxed_status, " relaxed_LB   ",relaxed_LB, "  UB  ",UB)
                    if relaxed_status == :Optimal && relaxed_LB >= UB
                            if Rprobing.redCosts[i] <= -1e-4
                                println("update P.colLower")
                                P.colLower[i] = xl + (relaxed_LB-UB)/abs(Rprobing.redCosts[i])
                                println("updated point", P.colLower[i] , "xl",xl, "   dual ", Rprobing.redCosts[i])
                            end
                    end
                    updateFirstBounds!(P, P.colLower, P.colUpper)

                    Rprobing.colLower[i] = xu
                    Rprobing.colUpper[i] = xu
                    Rprobing.colVal[i] = xu
                    relaxed_status = solve(Rprobing)
                    relaxed_LB = getRobjective(Rprobing, relaxed_status, LB)
                    println("relaxed_status  ", relaxed_status, " relaxed_LB   ",relaxed_LB, "  UB  ",UB)
                    if relaxed_status == :Optimal && relaxed_LB >= UB
                            if Rprobing.redCosts[i] <= -1e-4
                                println("update P.colUpper")
                                println("R    ", Rprobing.colVal[i])
                                P.colUpper[i] = xu - (relaxed_LB-UB)/abs(Rprobing.redCosts[i])
                                println("updated point", P.colUpper[i] , " xu ", xu, "   dual ", Rprobing.redCosts[i])
                            end
                    end
                    updateFirstBounds!(P, P.colLower, P.colUpper)





   		    xl = P.colLower[i]
    		    xu = P.colUpper[i]
		    if (xu - xl) <= small_bound_improve
		        i = i + 1
		    	continue    
		    end
		    updateStoBoundsFromSto!(P, Pprobing)		
		    leftPoint = xl + (xu - xl)/3
		    #=
		    if level == 1
		       leftPoint = xl + (xu - xl)/6
		    end
		    =#
		    updateFirstBounds!(Pprobing, xl, leftPoint, i)
		    updateExtensiveBoundsFromSto!(Pprobing, Pex_probing)
		    feasible_left = Sto_fast_feasibility_reduction!(Pprobing, pr_children, Pex_probing, prex, Rold, UB, LB, 0, true)	   
		    if !feasible_left
		        P.colLower[i] = leftPoint
		    else
			P.colLower[i] = Pprobing.colLower[i]	   
		    end	
		    updateFirstBounds!(P, P.colLower, P.colUpper)		    
		    updateStoBoundsFromSto!(P, Pprobing)


		    rightPoint = xu - (xu - xl)/3	
		    #=
		    if level ==	 1
		       rightPoint = xu - (xu - xl)/6
		    end
		    =#	
		    updateFirstBounds!(Pprobing, rightPoint, xu, i)
		    updateExtensiveBoundsFromSto!(Pprobing, Pex_probing)
		    feasible_right = Sto_fast_feasibility_reduction!(Pprobing, pr_children, Pex_probing, prex, Rold, UB, LB, 0, true)
		    if !feasible_right
		        P.colUpper[i] = rightPoint
		    else
			P.colUpper[i] = Pprobing.colUpper[i]	    
		    end
		    updateFirstBounds!(P, P.colLower, P.colUpper) 		    


		    #=
                    xl = P.colLower[i]
                    xu = P.colUpper[i]
                    updateExtensiveBoundsFromSto!(P, Pex_probing)
                    Rprobing = updaterelax(Rold, Pex_probing, prex, 1e10)
                    Rprobing.colCat[i] = :Fixed
                    Rprobing.colLower[i] = xl
                    Rprobing.colUpper[i] = xl
                    Rprobing.colVal[i] = xl
		    relaxed_status = solve(Rprobing)
                    relaxed_LB = getRobjective(Rprobing, relaxed_status, LB)
		    println("relaxed_status  ", relaxed_status, " relaxed_LB   ",relaxed_LB, "  UB  ",UB)
                    if relaxed_status == :Optimal && relaxed_LB >= UB             
                            if Rprobing.redCosts[i] <= -1e-4
                                println("update P.colLower")
                                P.colLower[i] = xl + (relaxed_LB-UB)/abs(Rprobing.redCosts[i])
                                println("updated point", P.colLower[i] , "xl",xl, "   dual ", Rprobing.redCosts[i])
                            end
                    end
                    updateFirstBounds!(P, P.colLower, P.colUpper)

		    Rprobing.colLower[i] = xu
                    Rprobing.colUpper[i] = xu
                    Rprobing.colVal[i] = xu		    
                    relaxed_status = solve(Rprobing)
                    relaxed_LB = getRobjective(Rprobing, relaxed_status, LB)
		    println("relaxed_status  ", relaxed_status, " relaxed_LB   ",relaxed_LB, "  UB  ",UB)
                    if relaxed_status == :Optimal && relaxed_LB >= UB
                            if Rprobing.redCosts[i] <= -1e-4
                                println("update P.colUpper")
                                println("R    ", Rprobing.colVal[i])
                                P.colUpper[i] = xu - (relaxed_LB-UB)/abs(Rprobing.redCosts[i])
                                println("updated point", P.colUpper[i] , " xu ", xu, "   dual ", Rprobing.redCosts[i])
                            end
                    end
                    updateFirstBounds!(P, P.colLower, P.colUpper)
		    =#


		    if !feasible_left || !feasible_right
			   n_reduced_probing += 1		      			   
                    	   println("A ", P.colLower, "     ", P.colUpper)
                    	   feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, LB, 0, true)
                    	   println("B ", P.colLower, "     ", P.colUpper, feasible)
			   if !feasible
			       n_reduced_probing = 0 
			       break
			   end			   			  
		    else
		        i = i + 1
		    end		    
	        end	
		#end
            	updateExtensiveBoundsFromSto!(P, Pex)
            	println("after feasibility reduction ", P.colLower, "   ", P.colUpper)
            	println("probing time:   ",  toq(), " (s)", feasible)
	    end

            if !feasible
               node_LB = UB 
               break
            end

	    # Lowerbound	    
	    tic()
	    relaxed_status, relaxed_LB, R, reduction_relax = getRelaxLowerBoundBT!(P, pr_children, Rold, Pex, prex, UB, node.LB)
	    println("Relax Lower time:   ",  toq(), " (s)")
            println("reduction_relax", reduction_relax)
            node_LB = max(node_LB,relaxed_LB)

	    if relaxed_status == :Optimal
	        #Rold.colVal = copy(R.colVal)
		if  ((UB-relaxed_LB)<= mingap || (UB-relaxed_LB) <= mingap*abs(relaxed_LB))
                    if ((UB-relaxed_LB)>=0) && (relaxed_LB <= FLB)
                        FLB = relaxed_LB
                    end
                    relaxed_status = :Infeasible
		end    
            end
            if relaxed_status == :Infeasible
	        node_LB = UB
                break
            end

	    node.LB = node_LB
            LB, ~ = getGlobalLowerBound(nodeList, FLB, UB)
            if LB >= node_LB
                LB = node_LB
            end

	end    

	    if (!feasible)  || (relaxed_status == :Infeasible)
                node_LB = UB
                continue
            end
	    
	    Rsol = copy(R.colVal)
	    #To do: check if can be moved within loop
	    if relaxed_status == :Optimal
	        Rold = copyModel(R)
	    end   
	    
	    WS_LB = -1e8
            WS_status = :NotOptimal
            WSfirst = zeros(nfirst)
            WSSol = nothing	    
	    tic()
            updateStoBoundsFromExtensive!(Pex, PWS)
            WS_status, WS_LB, WSfirst, WSSol = getWSLowerBound(PWS, node.parWSSol, UB, (UB-LB)/2, (UB-LB)/abs(LB)/2)
	    updateStoBoundsFromSto!(PWS, P) 
	    updateExtensiveBoundsFromSto!(P, Pex)
	    println("WS_status  ", WS_status, "  WS_LB ", WS_LB)
	    println("WS Lower time:   ",  toq(), " (s)")
            if WS_status == :Optimal && ( (UB-WS_LB)<= mingap || (UB-WS_LB)/abs(WS_LB) <= mingap)
                if ((UB-WS_LB)>=0) && (WS_LB <= FLB)
                    FLB = WS_LB
                end
                WS_status == :Infeasible
            end
            if WS_status == :Infeasible
	        node_LB = UB
		LB_UB_list[iter][3] = 2
                continue   #break
            elseif WS_status == :Optimal
                for idx in 1:nscen
                    v_map = Pex.ext[:v_map][idx]
                    scenarios[idx].colVal = copy(scenariosWS[idx].colVal)
                    scenarios[idx].colLower[end] = max(scenarios[idx].colLower[end], WSSol.secondobjVals[idx])
                    Pex.colLower[v_map[end]] = scenarios[idx].colLower[end]
                end
                P.colVal = copy(WSfirst)
            end
	   

	    LB_status = (WS_status == :Optimal || relaxed_status == :Optimal)? (:Optimal): (:NotOptimal)
	    node_LB_updated = false
	    if (WS_LB - node_LB) >= machine_error
	        node_LB = WS_LB
		node_LB_updated = true
		LB_UB_list[iter][3] = 2 
	    end
	    if LB_status == :Optimal		
		updateVarSelector(vs, node, P, node_LB)
	    end	
	    node.LB = node_LB
	    LB, ~ = getGlobalLowerBound(nodeList, FLB, UB)    
            if LB >= node_LB
                LB = node_LB
            end
	    println("FLB: ", FLB)
	    println("node L:", node_LB)
	    println("global L:", LB)
            #if lower bound is inferior than Upper bound, delete
            if ((UB-node_LB)<= mingap || ((UB-node_LB)/abs(node_LB)<= mingap))
	       continue  #break
            end


	    #Upper Bound
	    UB_status = :NotOptimal
	    UB_updated = false
	    tic()
	    UB_status, WSfix_status, node_UB, local_UB = getUpperBound!(Pex, prex, PWSfix, P, node_LB, WS_status, WSfirst, (UB-LB)/2, (UB-LB)/abs(LB)/2, level)
	    if (local_UB - node_UB) >= machine_error
	       LB_UB_list[iter][4] = 2
	    end

	    println("upper time:   ",  toq(), " (s)")
            #if better upperbound is achieved, update upper bound
            if (node_UB < UB)
	        UB_updated = true
      	        UB = node_UB
	    	x = copy(P.colVal)
	    	delete_nodes = []
	    	for (idx,n) in enumerate(nodeList)  	     
		    #println("new UB", idx,  "   ", node.LB)
	    	    if (((UB-n.LB)<= mingap) || ((UB-n.LB) <=mingap*abs(n.LB))) 
	                push!(delete_nodes, idx)
			if (((UB-n.LB)<= mingap) || ((UB-n.LB) <=mingap*abs(n.LB))) && ((UB-n.LB)>=0)
            		    if n.LB <= FLB
			        println("update FLB when better upperbound is achieved  ", FLB)
            		        FLB = n.LB
            		    end
	                end
	            end
                end
		deleteat!(nodeList, sort(delete_nodes))
            end	  
            println("UB:  ", UB)	   
	    
	    if UB_updated || node_LB_updated
	        feasible = Sto_medium_feasibility_reduction(P,pr_children, Pex, prex, Rold, UB, node_LB, bVarsId)
            	#feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, node_LB)
            	if !feasible
                    node_LB = UB
                    continue   #break
                end
            	updateExtensiveBoundsFromSto!(P, Pex)
            	println("after multiplier medium feasibility reduction ", P.colLower, "   ", P.colUpper)
	    end	    
    	#end

	#branch
	if ((UB-node_LB)<= mingap) || ((UB-node_LB) <= mingap*abs(node_LB))
	    if (UB-node_LB)>=0  && (node_LB <= FLB)
	        FLB = node_LB	     
		println("update FLB when closed to UB  ", FLB)
	    end
	else
	    	

	    max_score = 0
	    child_left_LB = node_LB
	    child_right_LB = node_LB	    
	    #=
	    if node.parmaxScore >= max_score_min   
	        tic()
		bVarId, max_score, exitBranch, child_left_LB, child_right_LB, FLB, node_LB = SelectVar!(vs, P, Rsol, WSfirst, WS_status, relaxed_status, PWS_child, Pex_child, WSSol, Rold, prex, UB, node_LB, Pex, FLB, node, pr_children, P_child)
		println("var select time:   ",  toq(), " (s)")
	    	if exitBranch
	            continue   
	        end
	    end		
	    if max_score <= max_score_min
                bVarId = SelectVarMaxRange(bVarsId, P) 
	    end
	    =#
	    bVarId = SelectVarMaxRange(bVarsId, P)				    
	    bValue = computeBvalue(Pex, P, bVarId, Rsol, WSfirst, WS_status, relaxed_status)
            println("bVarId:  ", bVarId, "  bValue: ", bValue, "        ", Pex.colLower[bVarId], "          ",Pex.colUpper[bVarId])
	    println("node_LB", node_LB)
            branch!(nodeList, P, Pex, bVarId, bValue, node, node_LB, WSSol, child_left_LB, child_right_LB, WS_status, relaxed_status, Rsol, max_score)

	    #=	    
            bVarId = SelectVarMaxRange(bVarsId, P)
            bValue = computeBvalue(Pex, P, bVarId, Rsol, WSfirst, WS_status, relaxed_status)
	    println("bVarId:  ", bVarId, "  bValue: ", bValue, "        ", Pex.colLower[bVarId], "          ",Pex.colUpper[bVarId])
	    branch!(nodeList, P, Pex, bVarId, bValue, node, node_LB, WSSol, -1e10, -1e10, WS_status, relaxed_status, Rsol, 0)	   
	    =#
	end
    end
    if nodeList==[]
       println("all node solved")
       if UB > 0
            LB = min(UB-mingap, UB/(1+mingap))
       else
	    LB = min(UB-mingap, UB/(1-mingap))
       end	    
    end
    if iter>=2 
        LB_UB_list[iter-1][1] = LB
        LB_UB_list[iter-1][2] = UB
    end   
    P.colVal = x
    P.objVal = UB
    println("Solution time:   ",  toq(), " (s)")
    println("solved nodes:  ",iter)
    @printf "%-52d  %-14.4e %-14.4e %-7.4f %s \n" iter  LB UB (UB-LB)/abs(LB)*100 "%"
    println("first stage sol   ",x)
    println("LB_UB_list   ", LB_UB_list)
    return P, LB_UB_list
end


function getRelaxLowerBoundBT!(P, pr_children, Rold, Pex, prex, UB, defaultLB, nsolve = 20)	 
	    reduction_all = 0
	    reduction_first = 0
	    nfirst = P.numCols
	    xlold = copy(Pex.colLower)
            xuold = copy(Pex.colUpper)	
            relaxed_status, relaxed_LB, R = getRelaxLowerBound(Rold, Pex, prex, UB, defaultLB, nsolve)
	    n_reduced_cost_BT = 0
            if relaxed_status == :Optimal && ( (UB-relaxed_LB)<= mingap || (UB-relaxed_LB)/abs(relaxed_LB) <= mingap)
	        return (relaxed_status, relaxed_LB, R, reduction_first)
            end
            if relaxed_status == :Infeasible
                relaxed_LB = UB
		return (relaxed_status, relaxed_LB, R, reduction_first)
            elseif relaxed_status == :Optimal
		fb = fast_feasibility_reduction!(R, nothing, UB)
            	if !fb
                   relaxed_status = :Infeasible
                   relaxed_LB = UB
		   return (relaxed_status, relaxed_LB, R, reduction_first)
                end
        	Pex.colLower = R.colLower[1:length(Pex.colLower)]
        	Pex.colUpper = R.colUpper[1:length(Pex.colUpper)]
                n_reduced_cost_BT = reduced_cost_BT!(Pex, prex, R, UB, relaxed_LB)
                updateStoBoundsFromExtensive!(Pex, P)
            end

	    left_all = 1
            for i in 1:length(Pex.colLower)
            	if (xuold[i] + Pex.colLower[i] - xlold[i] - Pex.colUpper[i]) > small_bound_improve
		   left_all = left_all * (Pex.colUpper[i] - Pex.colLower[i])/ (xuold[i]- xlold[i])
                end
       	    end
	    reduction_all = 1 - left_all

            if reduction_all > 0.1
	        bVarsId = collect(1:nfirst)
                #feasible = Sto_medium_feasibility_reduction(P, pr_children, Pex, prex, Rold, UB, relaxed_LB, bVarsId)
                feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, relaxed_LB)
                if !feasible
                   relaxed_status = :Infeasible
                   relaxed_LB = UB
                   return (relaxed_status, relaxed_LB, R, reduction_first)
                end
                updateExtensiveBoundsFromSto!(P, Pex)

            	left_first = 1
            	for i in 1:nfirst  #length(Pex.colLower)
                    if (xuold[i] + Pex.colLower[i] - xlold[i] - Pex.colUpper[i]) > small_bound_improve
                        left_first = left_first * (Pex.colUpper[i] - Pex.colLower[i])/ (xuold[i]- xlold[i])
                    end
                end
		reduction_first = 1 - left_first
                println("after multiplier medium feasibility reduction ", P.colLower, "   ", P.colUpper)
            end
	    return (relaxed_status, relaxed_LB, R, reduction_first)
end


# get lower bound from relaxed problem, and update bounds based on feasibility_reduction of R
function getRelaxLowerBound(Rold, Pex, prex, UB, defaultLB, nsolve =20)
            #solve relaxed problem and get lower bound of this node
            R = updaterelax(Rold, Pex, prex, UB)
	    R.solver = GurobiSolver(Method=1, Threads=1, LogToConsole=0,OutputFlag=0) 
            #R = relax(Pex, prex, UB)
            Rx = R.colVal
	    no_iter_wo_change = 0
            tic()
            relaxed_status = solve(R)
            relaxed_LB = getRobjective(R, relaxed_status, defaultLB)
            if relaxed_status == :Optimal
 	        Rx = R.colVal
                Ldecrease = 1e10
		#addOuterApproximationGrid!(R, prex, 20)
                for i=1:nsolve   #&& checkOuterApproximation(R, prex)  #while Ldecrease >= LP_improve_tol
                    newncon_convex = addOuterApproximation!(R, prex)
                    newncon_aBB = addaBB!(R, prex)
                    if (newncon_aBB + newncon_convex) != 0  
                        relaxed_status_trial = solve(R)
                        old_relaxed_LB = relaxed_LB
			relaxed_LB_trial = getRobjective(R, relaxed_status_trial, old_relaxed_LB)
			println("relaxed_LB_trial       ", relaxed_LB_trial)
                        if relaxed_status_trial == :Optimal
                            Rx = R.colVal
                            relaxed_LB = relaxed_LB_trial
                        elseif relaxed_status_trial == :Infeasible
                            relaxed_status = :Infeasible
                            relaxed_LB = UB
                            break
                        else
                            break
                        end
			Ldecrease = relaxed_LB_trial - old_relaxed_LB
			if Ldecrease >= LP_improve_tol
			    no_iter_wo_change = 0
			else
			    no_iter_wo_change += 1	    
			end
			if no_iter_wo_change >= 10
			    break
			end
                    else
                        break
                    end
                end
            end
            println("iterative relax time:   ",  toq(), " (s)", relaxed_status, "    ", relaxed_LB)
            R.colVal = Rx
            println("relaxed_LB  ",relaxed_LB)
            return (relaxed_status, relaxed_LB, R)
end


function getWSLowerBound(PWS::JuMP.Model, parWSSol, UB, absgap, regap)
            nfirst = PWS.numCols
            WSSol = StoSol(length(PlasmoOld.getchildren(PWS)))

            # get lower bound from wait and see formulation
            WS_LB = 0
            WS_status = :Optimal
            scenariosWS = PlasmoOld.getchildren(PWS)
            WSfirstSol = []
            for k = 1:nfirst
                push!(WSfirstSol, Float64[])
            end
            WSfirst = zeros(nfirst)
            nscen = length(scenariosWS)
            for (idx,scenario) in enumerate(scenariosWS)
                if parWSSol != nothing
                    secondSol = parWSSol.secondSols[idx]
                    scenario.colVal = copy(secondSol)
                    if inBounds(scenario, secondSol)
                        println("  in bounds  ")
                        for k = 1:nfirst
                            varid = scenario.ext[:firstVarsId][k]
                            if varid != -1
                               push!(WSfirstSol[k], scenario.colVal[varid])
                            end
                        end
                        WS_LB += parWSSol.secondobjVals[idx]     #getobjectivevalue(scenario)
                        WSSol.secondSols[idx] = secondSol
                        WSSol.secondobjVals[idx] =  parWSSol.secondobjVals[idx]
                        println("obj  ", parWSSol.secondobjVals[idx])
                        #println("sol  ", secondSol)
                        continue
                    end
                end

                scenariocopy = copyModel(scenario)
                scenariocopy.solver = IpoptSolver(print_level = 0)
                scenariocopy_status = solve(scenariocopy)
		# check unbounded, if unbounded, lowerbound->-1e10
		if (scenariocopy_status == :Unbounded) || (scenariocopy_status == :Optimal&& (getobjectivevalue(scenario) <= -1e10))
		    WS_status = :NotOptimal
                    WS_LB = - 1e20
                    break
                elseif scenariocopy_status == :Optimal && parWSSol == nothing ##
                    scenario.colVal = copy(scenariocopy.colVal)
                end
		
		regap = mingap
		absgap = mingap
                if scenariocopy_status == :Optimal
                    scenario.solver = SCIPSolver("limits/time", 100.0, "display/verblevel", 0, "limits/gap", max(mingap, regap)/2, "limits/absgap", max(mingap/nscen, absgap/nscen)/2, "setobjlimit", scenariocopy.objVal + 0.01)
                else
                    scenario.solver = SCIPSolver("limits/time", 100.0, "display/verblevel", 0, "limits/gap", max(mingap, regap)/2, "limits/absgap", max(mingap/nscen, absgap/nscen)/2)
                end
		#JuMP.build(scenario)
		#SCIP.setwarmstart!(scenario.internalModel, scenario.colVal)
                #scenario.solver = AmplNLSolver("/opt/scipoptsuite-3.2.1/scip-3.2.1/interfaces/ampl/bin/scipampl")
                #scenario.solver = BaronSolver()
                #push!(scenario.solver.options, (:EpsA, mingap/nscen), (:EpsR, mingap), (:MaxTime, 100))#,(:PrLevel, 0))

                scenario_status = solve(scenario)
		scenario_LB = getobjbound(scenario) #getobjbound(scenario.internalModel)
		if scenariocopy_status == :Optimal && scenario_status == :Infeasible
		    scenario_status = :Optimal
		    scenario.colVal = copy(scenariocopy.colVal)
		    scenario.objVal = scenariocopy.objVal
		    scenario_LB = scenariocopy.objVal
                end
                projection!(scenario.colVal, scenario.colLower, scenario.colUpper)
                println("  scenario :  ", idx, "   WS  scip status", scenario_status, "  obj ", getobjectivevalue(scenario), " lb  ", scenario_LB)
                #println("WS  scip solution     ", scenario.colVal)


                if scenario_status == :Optimal || scenario_status == :UserLimit
                    for k = 1:nfirst
                        varid = scenario.ext[:firstVarsId][k]
                        if varid != -1
                            push!(WSfirstSol[k], scenario.colVal[varid])
                        end
                    end
		    JuMP.setlowerbound(scenario[:objective_value], scenario_LB)
		    #=
                    if scenario_status == :UserLimit
                        println(scenario)
			error("WS_status UserLimit")
                    end
		    =#
                    WS_LB += scenario_LB                        #getobjectivevalue(scenario)
                    WSSol.secondSols[idx] = scenario.colVal
                    WSSol.secondobjVals[idx] = scenario_LB

		    
		    currentLB = 0
		    for (idx,scenariotemp) in enumerate(scenariosWS)
		        currentLB += JuMP.getlowerbound(scenariotemp[:objective_value])
		    end
		    if ( (UB-currentLB)<= mingap || (UB-currentLB) <= mingap*abs(currentLB))
		        WS_status = :Infeasible
		        break
                    end
                elseif   scenario_status == :Infeasible
                    WS_status = :Infeasible
                    break
                else
                    WS_status = :NotOptimal
                    println("lower bound scip not optimal")
		    error("WS_status  NotOptimal")
                    #break
                end
            end
	    if  WS_status == :Optimal
            	for k = 1:nfirst
                    WSfirst[k] = median(WSfirstSol[k])
                end
	    end
            return (WS_status, WS_LB, WSfirst, WSSol)
end



function getLowerBound(PWS::JuMP.Model, parWSSol, Rold, Pex, prex, UB, defaultLB)	    
	    nfirst = PWS.numCols
	    relaxed_status, relaxed_LB, R = getRelaxLowerBound(Rold, Pex, prex, UB, defaultLB)
	    if relaxed_status == :Optimal && ( (UB-relaxed_LB)<= mingap || (UB-relaxed_LB)/abs(relaxed_LB) <= mingap)
                return (:Optimal, relaxed_LB, :NotOptimal, :Optimal, relaxed_LB, nothing, R, nothing)
            end
	    
	    LB = relaxed_LB
	    LB_status = relaxed_status
            WS_LB = -1e10
            WS_status = :NotOptimal
            WSfirst = zeros(nfirst)
	    WSSol = nothing
	    updateStoBoundsFromExtensive!(Pex, PWS)
	    WS_status, WS_LB, WSfirst, WSSol = getWSLowerBound(PWS, parWSSol, UB, mingap, mingap)
	    
            if WS_status == :Infeasible
                return (:Infeasible, UB, 0, 0, 0, UB, 0, 0)
            end
            println("WS_status  ", WS_status, "  WS_LB ", WS_LB)

            LB_status = :NotOptimal
            if WS_status == :Optimal || relaxed_status == :Optimal
                LB_status = :Optimal
            end
            LB = max(relaxed_LB, WS_LB)
            return (LB_status, LB, WS_status, relaxed_status, relaxed_LB, WSfirst, R, WSSol)
end


function getUpperBound!(Pex, prex, PWSfix::JuMP.Model, P, node_LB, WS_status, WSfirst, absgap, regap, level)
	    updateStoBoundsFromSto!(P, PWSfix)				    
	    nfirst = P.numCols
            UB = 1e10
	    UB_status = false
	    # to do, check if the solution of relaxation and WS are feasible
 	    WSfix_status = :NotOptimal


            local_status = Ipopt_solve(P) #ParPipsNlp_solve(P)
	    local_UB = 1e10
            if local_status == :Optimal
                    local_UB = getsumobjectivevalue(P)
                    if local_UB < UB
                        UB_status = local_status
                        UB = local_UB
                    end
                    println("local_UB  ", local_status, local_UB)
            end

	    #=
	    local_status = :NotOptimal
	    updateExtensiveBoundsFromSto!(P, Pex)
            Pex.solver = IpoptSolver()
            local_status = solve(Pex)
	    local_UB = 1e10
            if local_status == :Optimal
	        updateStoSolFromExtensive!(Pex, P)
		local_UB = Pex.objVal
                if local_UB < UB
                    UB_status = local_status
                    UB = local_UB
                end
		println("local_UB  ", local_status, local_UB)
            end
	    =#

	    if level <=2 || level%3 == 0
		if ((local_status == :Optimal || WS_status==:Optimal) && (UB-node_LB)>= mingap && (UB-node_LB)/abs(node_LB) >=mingap)
		    if local_status == :Optimal
		        WSfixfirst = copy(P.colVal)  
		    else   
		    	WSfixfirst = copy(WSfirst)
		    end	
                    WSfix_UB = 0
                    WSfix_status = :Optimal
		    scenariosWSfix = PlasmoOld.getchildren(PWSfix)
		    scenarios = PlasmoOld.getchildren(P)
            	    nscen = length(scenariosWSfix)

                    for (idx,scenarioWSfix) in enumerate(scenariosWSfix)
                    	if local_status == :Optimal
                            scenarioWSfix.colVal = copy(scenarios[idx].colVal)
                        end
                        firstVarsId = scenarioWSfix.ext[:firstVarsId]
                   	for k = 1:nfirst
                            varid = firstVarsId[k]
                            if varid != -1
			        scenarioWSfix.colCat[varid] = :Fixed
                        	scenarioWSfix.colLower[varid] = WSfixfirst[k]
                        	scenarioWSfix.colUpper[varid] = WSfixfirst[k]
				scenarioWSfix.colVal[varid] = WSfixfirst[k]	
			    end
                        end			
						
			regap = mingap
			absgap = mingap
			if local_status == :Optimal			
			    scenarioWSfix.solver = SCIPSolver("display/verblevel", 0, "limits/gap", max(mingap, regap)/2, "limits/absgap", max(mingap/nscen,absgap/nscen)/2, "limits/time", 100.0, "setobjlimit", scenarios[idx].objVal + 0.01)
			else 
			    scenarioWSfix.solver = SCIPSolver("display/verblevel", 0, "limits/gap", max(mingap, regap)/2, "limits/absgap", max(mingap/nscen,absgap/nscen)/2, "limits/time")	     
			end
			#JuMP.build(scenarioWSfix)
			#SCIP.setwarmstart!(scenarioWSfix.internalModel, scenarioWSfix.colVal)
			#scenarioWSfix.solver = BaronSolver()
                        #push!(scenarioWSfix.solver.options, (:EpsA, mingap/nscen), (:EpsR, mingap), (:MaxTime, 100))#,(:PrLevel, 0))
                        scenarioWSfix_status = solve(scenarioWSfix)

			if local_status == :Optimal && scenarioWSfix_status == :Infeasible
			   scenarioWSfix_status = :Optimal
   			   scenarioWSfix.colVal = copy(scenarios[idx].colVal)
			   scenarioWSfix.objVal =  scenarios[idx].objVal			
                        end

			projection!(scenarioWSfix.colVal, scenarioWSfix.colLower, scenarioWSfix.colUpper)
                        #println(scenarioWSfix_status, "   obj     ", getobjectivevalue(scenarioWSfix), "  sol to WSfix  ", scenarioWSfix.colVal)
			println("scenario   ", idx, "    WSfix  scip status       ", scenarioWSfix_status, "  obj ", getobjectivevalue(scenarioWSfix), " lb  ", getobjbound(scenarioWSfix))
                	#println("WSfix  scip solution     ", scenarioWSfix.colVal)

                        if scenarioWSfix_status == :Optimal 
                            WSfix_UB += getobjectivevalue(scenarioWSfix)
                        elseif   scenarioWSfix_status == :Infeasible
                            WSfix_status = :Infeasible
                            break
                        else
                            WSfix_status = :NotOptimal
                            println("scip not optimal")
                            break
                        end
                    end
                    if WSfix_status == :Optimal 
			UB_status == WSfix_status
			if WSfix_UB < UB
                            UB = WSfix_UB
			end
                    end
                    println("WSfix_status  ", WSfix_status, "  WSfix_UB ", WSfix_UB)
                end	
	    end	
	    println("node U:", UB)
	    return UB_status, WSfix_status, UB, local_UB 
end


function SelectVar!(vs, P, Rsol, WSfirst, WS_status, relaxed_status, PWS_child, Pex_child, WSSol, Rold, prex, UB, node_LB, Pex, FLB, node, pr_children, P_child)
    node_LB_old = node_LB
    nfirst = length(PWS_child.colVal)
    exitBranch = false	 
    child_left_LB = -1e10
    child_right_LB = -1e10
    updateScore!(vs, Pex, P, Rsol, WSfirst, WS_status, relaxed_status)
    max_score = 0
    bVarId = vs.varId[1]
    no_consecutive_upates_wo_change = 0
    println("before select", P.colLower, P.colUpper)
    n_rel = 100
    lambda_consecutive = 8
    #for i in 1:length(vs.varId)
    i = 1
    ninfeasible = 0
    while i <= length(vs.varId)
            varId = vs.varId[i]
            bVarl = Pex.colLower[varId]
            bVaru = Pex.colUpper[varId]

	    if (bVaru - bVarl) <= small_bound_improve
	       i = i + 1
	       continue
	    end
	    bValue = computeBvalue(Pex, P, varId, Rsol, WSfirst, WS_status, relaxed_status)
            if min(vs.n_right[i], vs.n_left[i]) < n_rel && no_consecutive_upates_wo_change <= lambda_consecutive #&& (!vs.tried[i]) 
                vs.tried[i] = true
                updateStoBoundsFromExtensive!(Pex, PWS_child)
                Pex_child.colLower = copy(Pex.colLower)
                Pex_child.colUpper = copy(Pex.colUpper)
                updateFirstBounds!(PWS_child, bVarl, bValue, varId)
                updateExtensiveFirstBounds!(Pex_child, PWS_child, bVarl, bValue, varId)
		println("varId:   ",varId, "bValue:  ", bValue, "  xl  ", bVarl, "  xu  ", bVaru, "node_LB  ", node_LB)


		LB_status_left, obj_left, ~ = getRelaxLowerBound(Rold, Pex_child, prex, UB, node_LB)
                improve_left = 0
                if LB_status_left == :Optimal
                     if (((UB- obj_left)<= mingap) || ((UB- obj_left) <=mingap*abs(obj_left)))
                         if ((UB-obj_left)>=0) && (obj_left <= FLB)
                             FLB = obj_left
			     println("update FLB when select left var  ", FLB)
                         end
			 P.colLower[varId] = bValue
                         Pex.colLower[varId] = bValue
                         LB_status_left = :Infeasible
                     end
                    improve_left = max(0, obj_left - node_LB_old)
                    nlinfeasibility = bValue - bVarl
                    vs.pcost_left[i] = (vs.pcost_left[i]*vs.n_left[i] + improve_left/nlinfeasibility)/(vs.n_left[i] + 1)
                    vs.n_left[i] += 1
                elseif LB_status_left == :Infeasible
                    #if LB problem is infeasible, update bounds
                    P.colLower[varId] = bValue
		    Pex.colLower[varId] = bValue
		    improve_left = max(0, UB - node_LB_old)
                    #continue
                else
                    #if not solved sucessfully, and cannot be proven to be infeasible, continue
                    continue
                end
		println("lb  ",Pex.colLower[varId], "ub   ",Pex.colUpper[varId])
		println("LB_status_left:  ", LB_status_left, "  obj  ", obj_left)


		if LB_status_left == :Infeasible
		    updateStoBoundsFromExtensive!(Pex, P)
		    feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, node_LB)
		    updateExtensiveBoundsFromSto!(P, Pex)		    
		    updateExtensiveBoundsFromSto!(P, Pex_child)
		    updateStoBoundsFromSto!(P, PWS_child)
		end


                updateStoBoundsFromExtensive!(Pex, PWS_child)
                Pex_child.colLower = copy(Pex.colLower)
                Pex_child.colUpper = copy(Pex.colUpper)
                updateFirstBounds!(PWS_child, bValue, bVaru, varId)
                updateExtensiveFirstBounds!(Pex_child, PWS_child, bValue, bVaru, varId)

		LB_status_right, obj_right, ~ = getRelaxLowerBound(Rold, Pex_child, prex, UB, node_LB)
                improve_right = 0
		println("LB_status_right:  ", LB_status_right, "  obj  ", obj_right)
                if LB_status_right == :Optimal
                    if (((UB- obj_right)<= mingap) || ((UB- obj_right) <=mingap*abs(obj_right)))
                         if ((UB-obj_right)>=0) && (obj_right <= FLB)
                             FLB = obj_right
			     println("update FLB when select right  ", FLB)
                         end
                         Pex.colUpper[varId] = bValue
			 P.colUpper[varId] = bValue
                         LB_status_right = :Infeasible
			 if obj_right >= UB
			     obj_right = UB
			 end   
                    end
                    improve_right = max(0, obj_right - node_LB_old)
                    nlinfeasibility = bValue - bVarl
                    vs.pcost_right[i] = (vs.pcost_right[i]*vs.n_right[i] + improve_right/nlinfeasibility)/(vs.n_right[i] + 1)
                    vs.n_right[i] += 1
                    vs.score[i] = vs.pcost_right[i] * nlinfeasibility
                elseif LB_status_right == :Infeasible
                    #if relaxed problem is infeasible, update bounds
                    P.colUpper[varId] = bValue
                    Pex.colUpper[varId] = bValue
		    obj_right = UB
                    improve_right = max(0, obj_right - node_LB_old)
                    #continue
                else
                    #if not solved sucessfully, and cannot be proven to be infeasible, continue
                    continue
                end
		println("lb   ub",Pex.colLower[varId], Pex.colUpper[varId])
                println("LB_status_right:  ", LB_status_right, "  obj  ", obj_right)


                if LB_status_right == :Infeasible
		    updateStoBoundsFromExtensive!(Pex, P)
                    feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, node_LB)
                    updateExtensiveBoundsFromSto!(P, Pex)
                end


                if LB_status_right == :Infeasible && LB_status_left == :Infeasible
		        exitBranch = true
			break
                end		
		
		if (LB_status_right == :Infeasible || LB_status_left == :Infeasible)  #&& ninfeasible <= 1
		   ninfeasible += 1
		   continue
		end

		if LB_status_right == :Infeasible && LB_status_left == :Optimal && (obj_left-node_LB)>=0
		    node_LB = obj_left
		elseif  LB_status_left == :Infeasible && LB_status_right == :Optimal && (obj_right-node_LB)>=0
		    node_LB = obj_right	   
		elseif LB_status_left == :Optimal && LB_status_right == :Optimal && (obj_left-node_LB)>=0 && (obj_right-node_LB)>=0
		    node_LB = min(obj_left, obj_right)   
		end

                vs.score[i] = compute_score(improve_left, improve_right)
                println("improve_left:  ", improve_left, "   improve_right:  ", improve_right, " score ",vs.score[i])

		println("LB_status_left  ", LB_status_left, "  obj_left  ", obj_left, "  LB_status_right  ", obj_right)
        	if vs.score[i] >= max_score
                    max_score = vs.score[i]
                    bVarId = varId
                    if LB_status_left == :Optimal
                        child_left_LB = obj_left
                        if LB_status_right == :Infeasible
                    	    child_right_LB = child_left_LB
                        end
                    end
                    if LB_status_right == :Optimal
                        child_right_LB = obj_right
                        if LB_status_left == :Infeasible
                    	    child_left_LB = obj_right
                        end
                    end
                    no_consecutive_upates_wo_change = 0
                else
                    no_consecutive_upates_wo_change += 1
                end
            end
            if vs.score[i] >= max_score
                max_score = vs.score[i]
                bVarId = varId
            end
            if (no_consecutive_upates_wo_change > lambda_consecutive && max_score > 0)
                break
            end
	    i = i + 1
	    ninfeasible = 0
    end
    println("bVarId:  ", bVarId, "max_score  ", max_score)
    updateStoBoundsFromExtensive!(Pex, P)     
    println("after select", P.colLower, P.colUpper)
    return bVarId, max_score, exitBranch, child_left_LB, child_right_LB, FLB, node_LB
end


function branch!(nodeList, P, Pex, bVarId, bValue, node, node_LB, WSSol, child_left_LB, child_right_LB, WS_status, relaxed_status, Rsol, max_score)
    xl = copy(Pex.colLower)
    xu = copy(Pex.colUpper)
    println("node_LB in branch   ", node_LB)
    updateExtensiveFirstBounds!(xl, xu, P, Pex, Pex.colLower[bVarId], bValue, bVarId)
    if WS_status == :Optimal
        left_node = Node(xl, xu, bVarId, node.level+1, node_LB, -1, WSSol, nothing, max_score)
    else
        left_node = Node(xl, xu, bVarId, node.level+1, node_LB, -1, nothing, nothing, max_score)
    end
    if child_left_LB != -1e10 && (child_left_LB >= node_LB)
       left_node.LB = child_left_LB
       left_node.direction = 0
    end

    xl = copy(Pex.colLower)
    xu = copy(Pex.colUpper)
    updateExtensiveFirstBounds!(xl, xu, P, Pex, bValue, Pex.colUpper[bVarId], bVarId)
    if WS_status == :Optimal
        right_node = Node(xl, xu, bVarId, node.level+1, node_LB, 1, WSSol, nothing, max_score)
    else
        right_node = Node(xl, xu, bVarId, node.level+1, node_LB, 1, nothing, nothing,max_score)
    end
    if child_right_LB != -1e10 && (child_right_LB >= node_LB)
       right_node.LB = child_right_LB
       right_node.direction = 0
    end
    if relaxed_status == :Optimal
       left_node.parRSol = copy(Rsol)
       right_node.parRSol = copy(Rsol)
    end

    push!(nodeList, left_node)
    push!(nodeList, right_node)
    println("left ", nodeList[end-1].LB, "   ", left_node.LB)
    println("right ", nodeList[end].LB, "   ", right_node.LB) 
end


#=
typealias LinearConstraint GenericRangeConstraint{AffExpr}

type GenericRangeConstraint{TermsType} <: AbstractConstraint
    terms::TermsType
    lb::Float64
    ub::Float64
end
m.quadconstr
m.linconstr

# ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c  [≤,≥]  0
type GenericQuadConstraint{QuadType} <: AbstractConstraint
    terms::QuadType
    sense::Symbol
end
typealias QuadConstraint GenericQuadConstraint{QuadExpr}


type GenericQuadExpr{CoefType,VarType} <: AbstractJuMPScalar
    qvars1::Vector{VarType}
    qvars2::Vector{VarType}
    qcoeffs::Vector{CoefType}
    aff::GenericAffExpr{CoefType,VarType}
end

typealias QuadExpr GenericQuadExpr{Float64,Variable}
=#





