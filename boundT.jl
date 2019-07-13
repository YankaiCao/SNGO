function Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, U, L = -1e20, ngrid = 0, solve_relax = true)
    nfirst = P.numCols
    scenarios = PlasmoOld.getchildren(P)
    nscen = length(scenarios)
    feasibility_reduced = 1e10
    rn = 1
    changed = trues(nscen)
    updateStoFirstBounds!(P)

    while feasibility_reduced >= 1  
        xlold = copy(P.colLower)
        xuold = copy(P.colUpper)
	updateExtensiveBoundsFromSto!(P, Pex)

	#R = relax(Pex, prex, U)
	R = updaterelax(Rold, Pex, prex, U)
	if ngrid > 0
	    addOuterApproximationGrid!(R, prex, ngrid)
	end	
        fb = fast_feasibility_reduction!(R, nothing, U)
        if !fb
            return fb
        end
        Pex.colLower = R.colLower[1:length(Pex.colLower)]
        Pex.colUpper = R.colUpper[1:length(Pex.colUpper)]


	if solve_relax && rn == 1
	    relaxed_status = solve(R)
	    relaxed_LB = getRobjective(R, relaxed_status, L)
            if relaxed_status == :Optimal && ( (U-relaxed_LB)<= mingap || (U-relaxed_LB)/abs(relaxed_LB) <= mingap)
	        #FLB = relaxed_LB
	        return false
            elseif relaxed_status == :Infeasible
                return false
            end
            if relaxed_status == :Optimal
                n_reduced_cost_BT = reduced_cost_BT!(Pex, prex, R, U, relaxed_LB)
		#if n_reduced_cost_BT >= 1
                #    println("inside Sto_fast after reduced cost  ", Pex.colLower[1:nfirst], "   ", Pex.colUpper[1:nfirst])
		#end    
            end
	end
	updateStoBoundsFromExtensive!(Pex, P)
	updateStoFirstBounds!(P)


        for (idx,scenario) in enumerate(scenarios)
	    #println("scenario",idx)
	    if changed[idx]
                fb = fast_feasibility_reduction!(scenario, pr_children[idx], 1e10)
            	if !fb
                    return fb
                end
	    end	
        end


        for (idx,scenario) in enumerate(scenarios)
            firstVarsId = scenario.ext[:firstVarsId]
            for i in 1:nfirst
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

	for (idx,scenario) in enumerate(scenarios)
	    changed[idx] = false
            for j = 1:nfirst
	    	varid = scenario.ext[:firstVarsId][j]
                if varid != -1
#		    if (P.colLower[j] - scenario.colLower[varid]) >= machine_error || (scenario.colUpper[varid]-P.colUpper[j]) >= machine_error
			scenario.colLower[varid] = P.colLower[j]
                    	scenario.colUpper[varid] = P.colUpper[j]
		    	changed[idx] = true
#		    end	
                end
            end
        end
        #updateFirstBounds!(P, P.colLower, P.colUpper)

	xls = Array{Float64}(nscen)
	xus = Array{Float64}(nscen)
        for j in 1:nscen
 	    xls[j] = getlowerbound(scenarios[j][:objective_value])
	    xus[j] = getupperbound(scenarios[j][:objective_value])
        end
	linearBackward!(ones(nscen), xls, xus, sum(xls), sum(xus),  L,  U )
	
        feasibility_reduced = 0
        for j in 1:nscen
            var = scenarios[j][:objective_value]
            xu = getupperbound(var)
            xu_trial = xus[j]
            if (xu-xu_trial) >= probing_improve
                feasibility_reduced += 1
                setupperbound(var, xu_trial)
                changed[j] = true
                #println("obj:  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
            end
            xl = getlowerbound(var)
            xl_trial = xls[j]
            if (xl_trial-xl) >= probing_improve
                feasibility_reduced += 1
                setlowerbound(var, xl_trial)
                changed[j] = true
            end
        end
	#=
	feasibility_reduced = 0
    	for j in 1:nscen
            sum_min = 0
            for k in 1:nscen
                if k != j	      
                    xlk = getlowerbound(scenarios[k][:objective_value])
                    sum_min += xlk
                end
            end
            var = scenarios[j][:objective_value]
            xu = getupperbound(var)
            xu_trial = (U - sum_min)
            if (xu-xu_trial) >= machine_error
                feasibility_reduced += 1
                setupperbound(var, xu_trial)
		changed[j] = true
                #println("obj:  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
            end
        end
	=#

        for i in 1:nfirst
	    if (xuold[i] + P.colLower[i] - xlold[i] - P.colUpper[i]) > small_bound_improve
                feasibility_reduced = feasibility_reduced + 1
            end
        end

        #println("feasibiltiy_based reduction round: ",rn, "changed bound: ", feasibility_reduced)
        rn += 1
    end
    #println("bye")
    return true
end





function fast_feasibility_reduction!(P, pr, U)
    n = P.numCols
    feasible = true
    feasibility_reduced = 1e10
    rn = 1
    if sum(P.colLower.<=P.colUpper) < n
       feasible = false
       return feasible
    end
    
    while feasibility_reduced >= 1  #0.01*ncols
       	xlold = copy(P.colLower)
	xuold = copy(P.colUpper)
	feasible = fast_feasibility_reduction_inner!(P, pr, U)
	if feasible == false
	    break
	end
	feasibility_reduced = 0
	for i in 1:n
	    if P.colLower[i] - xlold[i] >= small_bound_improve
	     	feasibility_reduced = feasibility_reduced + 1
	    end
	    if xuold[i] - P.colUpper[i] >= small_bound_improve
                feasibility_reduced = feasibility_reduced + 1
	    end	 
         end
         #println("feasibiltiy_based reduction round: ",rn, "changed bound: ", feasibility_reduced)
         rn += 1
    end
    return feasible
end



function linearBackward!(coeffs, xls, xus, sum_min, sum_max,  lb,  ub )
	  changed = false
          if sum_min >= lb && sum_max <=ub
              return changed
          end
    	  n_reduced = 1e10
    	  while n_reduced >= 1
	      n_reduced = 0
              for j in 1:length(xls)
                  alpha = coeffs[j]
              	  xl = xls[j]
              	  xu = xus[j]
              	  if alpha >= 0
                      sum_min_except = sum_min - alpha*xl
                      sum_max_except = sum_max - alpha*xu
              	  else
		      sum_min_except = sum_min - alpha*xu
                      sum_max_except = sum_max - alpha*xl
                  end
              	  xu_trial = xu
              	  xl_trial = xl
              	  if alpha > 0
                      xu_trial = (ub - sum_min_except)/alpha
                      xl_trial = (lb - sum_max_except)/alpha
                  elseif alpha < 0
                      xl_trial = (ub - sum_min_except)/alpha
                      xu_trial = (lb - sum_max_except)/alpha
              	  end
		  
              	  if (xl_trial-xl) >= (machine_error)
                      xls[j] = xl_trial
		      changed = true
		      n_reduced += 1
                      #println("col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
                  end
              	  if (xu-xu_trial) >= (machine_error)
                      xus[j] = xu_trial
		      changed = true
		      n_reduced += 1
                      #println("col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
              	  end
              end
	  end
	  return changed
end


function AffineBackward!(aff, sum_min, sum_max,  lb,  ub )
    xls = getlowerbound(aff.vars)
    xus = getupperbound(aff.vars)
    changed  = linearBackward!(copy(aff.coeffs), xls, xus, sum_min, sum_max,  lb,  ub)
    #if changed
        for j in 1:length(aff.vars)
            var = aff.vars[j]
            setlowerbound(var, xls[j])
            setupperbound(var, xus[j])
        end
    #end	  
end

function multiVariableForward(mv::MultiVariable)
    terms = mv.terms	
    qvars1 = terms.qvars1
    qvars2 = terms.qvars2
    qcoeffs = terms.qcoeffs
    aff = terms.aff
    
    #println("terms   ", terms)
    sum_min, sum_max = Interval_cal(terms)
    #println(sum_min, "     ", sum_max)



    varsincon = [aff.vars;qvars1;qvars2]
    varsincon = union(varsincon)
    # square_coeff x^2 + [alpha_min, alpha_max]x + [sum_min,sum_max]
    for var in varsincon
              #println("var  ", var.col, "  ",getlowerbound(var), "      ",getupperbound(var))
              square_coeff = 0
              alpha_min = 0
              alpha_max = 0
              sum_min_trial = aff.constant
              sum_max_trial = aff.constant

              for k in 1:length(aff.vars)
                  if var == aff.vars[k]
                     alpha_min += aff.coeffs[k]
                     alpha_max += aff.coeffs[k]
                  else
                     xlk = getlowerbound(aff.vars[k])
                     xuk = getupperbound(aff.vars[k])
                     coeff = aff.coeffs[k]
                     sum_min_trial += min(coeff*xlk, coeff*xuk)
                     sum_max_trial += max(coeff*xlk, coeff*xuk)
                  end
              end
              for k in 1:length(qvars1)
                  if qvars1[k] == var && qvars2[k] == var
                      square_coeff += qcoeffs[k]
                  elseif qvars1[k] != var && qvars2[k] != var
                      xlk1 = getlowerbound(qvars1[k])
                      xuk1 = getupperbound(qvars1[k])
                      xlk2 = getlowerbound(qvars2[k])
                      xuk2 = getupperbound(qvars2[k])
                      coeff = qcoeffs[k]
                      if qvars1[k] == qvars2[k]
                         if xlk1 <= 0 && 0 <= xuk1
                            sum_min_trial += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
                            sum_max_trial += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
                         else
                            sum_min_trial += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                            sum_max_trial += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                         end
                      else
                         sum_min_trial += min(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
                         sum_max_trial += max(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
                      end
                elseif qvars1[k] == var
                      xlk2 = getlowerbound(qvars2[k])
                      xuk2 = getupperbound(qvars2[k])
                      coeff = qcoeffs[k]
                      alpha_min += min(coeff*xlk2, coeff*xuk2)
                      alpha_max += max(coeff*xlk2, coeff*xuk2)
                  else
                      xlk1 = getlowerbound(qvars1[k])
                      xuk1 = getupperbound(qvars1[k])
                      coeff = qcoeffs[k]
                      alpha_min += min(coeff*xlk1, coeff*xuk1)
                      alpha_max += max(coeff*xlk1, coeff*xuk1)
                  end
             end
	     # square_coeff x^2 + [alpha_min, alpha_max]x + [sum_min,sum_max]

             xl = getlowerbound(var)
             xu = getupperbound(var)

	     if abs(square_coeff) == 0
                sum_min_temp = sum_min_trial + min(alpha_min*xl, alpha_min*xu, alpha_max*xl, alpha_max*xu)
                sum_max_temp = sum_max_trial + max(alpha_min*xl, alpha_min*xu, alpha_max*xl, alpha_max*xu)

             	if sum_min_temp > sum_min
                    sum_min = sum_min_temp
                end
             	if sum_max_temp < sum_max
                    sum_max = sum_max_temp
                end

            end
	    #=
	    if abs(square_coeff) > 0
                sum_min_temp = sum_min_trial
                sum_max_temp = sum_max_trial
                if alpha_min <= 0 && alpha_max >=0
                   sum_min_temp += min(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4, 0)
                   sum_max_temp += max(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4, 0)
                else
                   sum_min_temp += min(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4)
                   sum_max_temp += max(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4)
                end
		add sq*(x-b/2sq)^2	
	     	if sum_min_temp > sum_min
	     	   sum_min = sum_min_temp
	        end
	     	if sum_max_temp < sum_max
	     	   sum_max = sum_max_temp
	        end	
             end 
	     =#
    end	     



    # a x^2 + b x
    if (length(qvars1) == 1) 
        if qvars1[1] == qvars2[1]
	    a = qcoeffs[1]
	    var = qvars1[1]
	    @assert a != 0
	    @assert aff.constant == 0
	    @assert length(aff.vars) <= 1
	    b = 0
	    if length(aff.coeffs) == 1
	        b = aff.coeffs[1]
	    end   
            xl = getlowerbound(var)
            xu = getupperbound(var)
            sum_min = min(a*xl^2+ b*xl, a*xu^2+ b*xu)
            sum_max = max(a*xl^2+ b*xl, a*xu^2+ b*xu)
            if xl <= (-b/2a) && (-b/2a) <= xu
                sum_min = min(sum_min, -b^2/a/4.0)
		sum_max = max(sum_max, -b^2/a/4.0)
            end
	end
	#println(sum_min, "    ",sum_max)
    end
 	 
    return (sum_min, sum_max)
end


function multiVariableBackward!(mv::MultiVariable, sum_min, sum_max, lb, ub)
    if sum_min >= lb && sum_max <=ub
        return
    end
    terms = mv.terms
    qvars1 = terms.qvars1
    qvars2 = terms.qvars2
    qcoeffs = terms.qcoeffs
    aff = terms.aff

    #println("backward!    ", sum_min,"   ", sum_max, "  ",lb, "  ",ub)
    #println(terms)

    varsincon = [aff.vars;qvars1;qvars2]
    varsincon = union(varsincon)
    # square_coeff x^2 + [alpha_min, alpha_max]x + [sum_min,sum_max]
    for var in varsincon

    	      #println("var  ", var.col, "  ",getlowerbound(var), "      ",getupperbound(var))
              square_coeff = 0
              alpha_min = 0
              alpha_max = 0
              sum_min = aff.constant
              sum_max = aff.constant

              for k in 1:length(aff.vars)
                  if var == aff.vars[k]
                     alpha_min += aff.coeffs[k]
                     alpha_max += aff.coeffs[k]
                  else
                     xlk = getlowerbound(aff.vars[k])
                     xuk = getupperbound(aff.vars[k])
                     coeff = aff.coeffs[k]
                     sum_min += min(coeff*xlk, coeff*xuk)
                     sum_max += max(coeff*xlk, coeff*xuk)
                  end
              end

              for k in 1:length(qvars1)
                  if qvars1[k] == var && qvars2[k] == var
                      square_coeff += qcoeffs[k]
                  elseif qvars1[k] != var && qvars2[k] != var
                      xlk1 = getlowerbound(qvars1[k])
                      xuk1 = getupperbound(qvars1[k])
                      xlk2 = getlowerbound(qvars2[k])
                      xuk2 = getupperbound(qvars2[k])
                      coeff = qcoeffs[k]
                      if qvars1[k] == qvars2[k]
                         if xlk1 <= 0 && 0 <= xuk1
                            sum_min += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
                            sum_max += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
                         else
                            sum_min += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                            sum_max += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                         end
                      else
                         sum_min += min(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
                         sum_max += max(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
                      end
                  elseif qvars1[k] == var
                      xlk2 = getlowerbound(qvars2[k])
                      xuk2 = getupperbound(qvars2[k])
                      coeff = qcoeffs[k]
                      alpha_min += min(coeff*xlk2, coeff*xuk2)
                      alpha_max += max(coeff*xlk2, coeff*xuk2)
                  else
                      xlk1 = getlowerbound(qvars1[k])
                      xuk1 = getupperbound(qvars1[k])
                      coeff = qcoeffs[k]
                      alpha_min += min(coeff*xlk1, coeff*xuk1)
                      alpha_max += max(coeff*xlk1, coeff*xuk1)
                  end
             end

             xl = getlowerbound(var)
             xu = getupperbound(var)
             sum_min_temp = sum_min
             sum_max_temp = sum_max
             if abs(square_coeff) > 0
                if (xl<=0)&& (xu>=0)
                   sum_min_temp += min(square_coeff*xl*xl, square_coeff*xu*xu, 0)
                   sum_max_temp += max(square_coeff*xl*xl, square_coeff*xu*xu, 0)
                else
                   sum_min_temp += min(square_coeff*xl*xl, square_coeff*xu*xu)
                   sum_max_temp += max(square_coeff*xl*xl, square_coeff*xu*xu)
                end
             end

             #println("square_coeff   ", square_coeff)
             #println("alpha min: ", alpha_min, " max: ", alpha_max)
             #println("sum  min: ", sum_min, "  sum_max: ", sum_max)
             #println("sum_temp  min: ", sum_min_temp, "  sum_max: ", sum_max_temp)

             xu_trial = xu
             xl_trial = xl

             if alpha_min > 0
                xu_trial = max((ub - sum_min_temp)/alpha_min, (ub - sum_min_temp)/alpha_max)
                xl_trial = min((lb - sum_max_temp)/alpha_min, (lb - sum_max_temp)/alpha_max)
             elseif alpha_max < 0
                xu_trial = max((lb - sum_max_temp)/alpha_min, (lb - sum_max_temp)/alpha_max)
                xl_trial = min((ub - sum_min_temp)/alpha_min, (ub - sum_min_temp)/alpha_max)
             end
             #println("trial:  ", xl_trial,"  ",xu_trial)
             xu_trial = min(xu_trial, xu)
             xl_trial = max(xl_trial, xl)


             #second round if sq != 0
             if abs(square_coeff) > 0
                sum_min_temp = sum_min
                sum_max_temp = sum_max
                sum_min_temp += min(alpha_min*xl_trial, alpha_min*xu_trial, alpha_max*xl_trial, alpha_max*xu_trial)
                sum_max_temp += max(alpha_min*xl_trial, alpha_min*xu_trial, alpha_max*xl_trial, alpha_max*xu_trial)
                sqrt_ub = sqrt(max((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff, 0))
			sqrt_lb = sqrt(max(min((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff),0))
                #println("suqare_lb", sqrt_lb, "  ", sqrt_ub)
                xu_trial = min(xu_trial, sqrt_ub)
                xl_trial = max(xl_trial, -sqrt_ub)
                if xl_trial<= -sqrt_lb && xu_trial<= sqrt_lb && xu_trial>= -sqrt_lb
                      xu_trial = - sqrt_lb
                elseif xu_trial >= sqrt_lb && xl_trial<= sqrt_lb && xl_trial>= -sqrt_lb
                      xl_trial = sqrt_lb
                end
                #println("square_trial:  ", xl_trial,"  ",xu_trial)
             end

             #third round if sq != 0
             if abs(square_coeff) > 0
	     	sum_min_temp = sum_min
                sum_max_temp = sum_max

                #println("sum2 ", sum_min_temp, "  ",sum_max_temp)
                #println("alpha_min", alpha_min)
		#println("alpha_max", alpha_max)
                #println("square_coeff", square_coeff)
                if alpha_min <= 0 && alpha_max >=0
                   sum_min_temp += min(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4, 0)
                   sum_max_temp += max(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4, 0)
                else
                   sum_min_temp += min(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4)
                   sum_max_temp += max(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4)
                end

                #println("sum", sum_min_temp, "  ",sum_max_temp)
                sqrt_ub = sqrt(max((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff, 0))
                sqrt_lb = sqrt(max(min((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff),0))
                #println("suqare_lb", sqrt_lb, "  ", sqrt_ub)

		temp_min = min(alpha_min/square_coeff/2, alpha_max/square_coeff/2)
                temp_max = max(alpha_min/square_coeff/2, alpha_max/square_coeff/2)

                xu_trial = min(xu_trial, sqrt_ub - temp_min)
                xl_trial = max(xl_trial, - sqrt_ub - temp_max)
                #x+temp_max > sqrt_lb or x+temp_min <=-sqrt_lb
                if (sqrt_lb-temp_max) > (-sqrt_lb - temp_min)
                   if xl_trial<= (-sqrt_lb - temp_min)  && xu_trial<= (sqrt_lb-temp_max) && xu_trial>= (-sqrt_lb - temp_min)
		            xu_trial = (-sqrt_lb - temp_min)
                   elseif xu_trial >= (sqrt_lb-temp_max)  && xl_trial<= (sqrt_lb-temp_max)   && xl_trial>= (-sqrt_lb - temp_min)
                      xl_trial = (sqrt_lb-temp_max)
                   end
                end
                #println("square_trial2:  ", xl_trial,"  ",xu_trial)
             end


             if (xl_trial-xl) >= machine_error
                  setlowerbound(var, xl_trial)
                  #println("quard   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
             end

	     if (xu-xu_trial) >= machine_error
                  setupperbound(var, xu_trial)
                  #println("quard  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
             end
    end


    
    # a x^2 + b x
    if (length(qvars1) == 1)
        if qvars1[1] == qvars2[1]
            a = qcoeffs[1]
            var = qvars1[1]
            @assert a != 0
            @assert aff.constant == 0
            @assert length(aff.vars) <= 1
            b = 0
            if length(aff.coeffs) == 1
                b = aff.coeffs[1]
            end
            xl = getlowerbound(var)
            xu = getupperbound(var)

            sqrt_ub = sqrt(max((ub + b^2/a/4.0)/a, (lb + b^2/a/4.0)/a, 0))
            sqrt_lb = sqrt(max(min((ub + b^2/a/4.0)/a, (lb + b^2/a/4.0)/a),0))

            temp = b/a/2
            xu_trial = min(xu, sqrt_ub - temp)
            xl_trial = max(xl, - sqrt_ub - temp)

            if sqrt_lb > 0
                   if xl_trial<= (-sqrt_lb - temp)  && xu_trial<= (sqrt_lb-temp) && xu_trial>= (-sqrt_lb - temp)
                            xu_trial = (-sqrt_lb - temp)
                   elseif xu_trial >= (sqrt_lb-temp)  && xl_trial<= (sqrt_lb-temp)   && xl_trial>= (-sqrt_lb - temp)
                      xl_trial = (sqrt_lb-temp)
                   end
            end

            if (xl_trial-xl) >= machine_error
                  setlowerbound(var, xl_trial)
                  #println("quard   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
            end

            if (xu-xu_trial) >= machine_error
                  setupperbound(var, xu_trial)
                  #println("quard  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
            end
        end
    end
    
end


function fast_feasibility_reduction_inner!(P, pr, U)
    #println(P.colLower)
    #println(P.colUpper)
    feasible = true 	 
    for i = 1:length(P.linconstr)
      	  con = P.linconstr[i] 
	  aff = con.terms
          lb = con.lb
	  ub = con.ub

	  # check if constraint is feasible
	  sum_min, sum_max = Interval_cal(aff)
	  if (sum_min-ub) >= machine_error || (lb-sum_max) >= machine_error
	     feasible = false	
	     #println(con)
	     #println("infeasible linear   ", i, "   ", sum_min, "  ", sum_max, "  ",lb,"   ",ub)     
	     #println(P.colLower)
	     #println(P.colUpper)
	     return (feasible)
	  end         
	  AffineBackward!(aff, sum_min, sum_max,  lb,  ub )
    end
    
    if length(P.quadconstr) > 0
    multiVariable_list = pr.multiVariable_list 
    for i = 1:length(P.quadconstr)
          con = P.quadconstr[i]
          terms = con.terms
          qvars1 = terms.qvars1
          qvars2 = terms.qvars2
          qcoeffs = terms.qcoeffs
          aff = terms.aff
	  #println(con)
	  #println(P.colLower)
    	  #println(P.colUpper)
          lb = 0
          ub = 0
	  if con.sense == :(<=)
	     lb = -1e20
	  end
	  if con.sense == :(>=)
	     ub = 1e20
	  end

	  #=
          sum_min, sum_max = Interval_cal(terms)
          if (sum_min-ub) >= machine_error || (lb-sum_max) >= machine_error
             #println(con)
             #println("infeasible quad   ", i, "   ", sum_min, "  ", sum_max, "  ",lb,"   ",ub)
             #println(P.colLower)
             #println(P.colUpper)
             feasible = false
             return (feasible)
          end
	  =#
	  #println(terms, "   ", lb, "   ",ub)	  	  


	  mv_con = copy(multiVariable_list[i], P)
	  #mv_con = multiVariable_list[i]
	  
	  mvs =	mv_con.mvs
	  remainaff = mv_con.aff
	  nmw = length(mv_con.mvs)
	  mv_sum_min = Array{Float64}(nmw + 1)
	  mv_sum_max = Array{Float64}(nmw + 1)	  
	  for j = 1:nmw
	      mv_sum_min[j], mv_sum_max[j] = multiVariableForward(mvs[j])
	  end
	  mv_sum_min[nmw+1], mv_sum_max[nmw+1] = Interval_cal(remainaff)
	  sum_min = sum(mv_sum_min)
	  sum_max = sum(mv_sum_max)
	  #println("remainaff:   ",remainaff)
	  #println(mv_sum_min[nmw+1], "    ", mv_sum_max[nmw+1] )
	  mv_lb = copy(mv_sum_min)
	  mv_ub = copy(mv_sum_max)

          changed = linearBackward!(ones(Float64, nmw+1), mv_lb, mv_ub, sum_min, sum_max,  lb,  ub)
          sum_min = sum(mv_lb)
          sum_max = sum(mv_ub)
	  #println("multi lb ub:         ", mv_lb, mv_ub)


	  # check if constraint is feasible
	  #sum_min, sum_max = Interval_cal(terms)
          if (sum_min-ub) >= machine_error || (lb-sum_max) >= machine_error
	     #println(con)
	     #println("infeasible quad   ", i, "   ", sum_min, "  ", sum_max, "  ",lb,"   ",ub)
             #println(P.colLower)
             #println(P.colUpper)
             feasible = false
             return (feasible)
          end
	  #println("mv_lb:      ", mv_lb, mv_ub)
	  #println(changed)

	  #if changed	  
          for j = 1:nmw
              multiVariableBackward!(mvs[j],mv_sum_min[j], mv_sum_max[j], mv_lb[j], mv_ub[j])
          end

	  AffineBackward!(remainaff, mv_sum_min[end], mv_sum_max[end],  mv_lb[end],  mv_ub[end])
          #println("after new:  ")
          #println(P.colLower)
          #println(P.colUpper)
	  

	  #=
	  # square_coeff x^2 + [alpha_min, alpha_max]x + [sum_min,sum_max]		
          varsincon = [aff.vars;qvars1;qvars2]
          varsincon = union(varsincon)
          for var in varsincon              
	      square_coeff = 0
	      alpha_min = 0
	      alpha_max = 0
              sum_min = aff.constant
              sum_max = aff.constant

              for k in 1:length(aff.vars)
                  if var == aff.vars[k]
		     alpha_min += aff.coeffs[k]
		     alpha_max += aff.coeffs[k]	       
		  else
                     xlk = getlowerbound(aff.vars[k])
                     xuk = getupperbound(aff.vars[k])
                     coeff = aff.coeffs[k]
                     sum_min += min(coeff*xlk, coeff*xuk)
                     sum_max += max(coeff*xlk, coeff*xuk)
                  end
              end

	      for k in 1:length(qvars1)
	      	  if qvars1[k] == var && qvars2[k] == var
		      square_coeff += qcoeffs[k]
		  elseif qvars1[k] != var && qvars2[k] != var  		     
		      xlk1 = getlowerbound(qvars1[k])
                      xuk1 = getupperbound(qvars1[k])
		      xlk2 = getlowerbound(qvars2[k])
                      xuk2 = getupperbound(qvars2[k])
                      coeff = qcoeffs[k]
		      if qvars1[k] == qvars2[k]
		         if xlk1 <= 0 && 0 <= xuk1 
			    sum_min += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
                            sum_max += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
			 else
			    sum_min += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                            sum_max += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                         end
		      else	
                      	 sum_min += min(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
                      	 sum_max += max(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)  	  
		      end   
		  elseif qvars1[k] == var    
		      xlk2 = getlowerbound(qvars2[k])
                      xuk2 = getupperbound(qvars2[k])
		      coeff = qcoeffs[k]
		      alpha_min += min(coeff*xlk2, coeff*xuk2)
		      alpha_max += max(coeff*xlk2, coeff*xuk2)
		  else
	              xlk1 = getlowerbound(qvars1[k])
                      xuk1 = getupperbound(qvars1[k])
		      coeff = qcoeffs[k]
		      alpha_min	+= min(coeff*xlk1, coeff*xuk1)
		      alpha_max += max(coeff*xlk1, coeff*xuk1)
		  end
	     end

             xl = getlowerbound(var)
             xu = getupperbound(var)
	     sum_min_temp = sum_min
	     sum_max_temp = sum_max
	     if abs(square_coeff) > 0
	     	if (xl<=0)&& (xu>=0)
		   sum_min_temp += min(square_coeff*xl*xl, square_coeff*xu*xu, 0)
		   sum_max_temp	+= max(square_coeff*xl*xl, square_coeff*xu*xu, 0) 	   
		else
		   sum_min_temp	+= min(square_coeff*xl*xl, square_coeff*xu*xu)
                   sum_max_temp += max(square_coeff*xl*xl, square_coeff*xu*xu)
		end
	     end

	     #println("square_coeff   ", square_coeff)
             #println("alpha min: ", alpha_min, " max: ", alpha_max)
	     #println("sum  min: ", sum_min, "  sum_max: ", sum_max)
	     #println("sum_temp  min: ", sum_min_temp, "  sum_max: ", sum_max_temp)

	     xu_trial = xu
	     xl_trial = xl

	     if alpha_min > 0
	     	xu_trial = max((ub - sum_min_temp)/alpha_min, (ub - sum_min_temp)/alpha_max)
                xl_trial = min((lb - sum_max_temp)/alpha_min, (lb - sum_max_temp)/alpha_max)	  
	     elseif alpha_max < 0
		xu_trial = max((lb - sum_max_temp)/alpha_min, (lb - sum_max_temp)/alpha_max) 
		xl_trial = min((ub - sum_min_temp)/alpha_min, (ub - sum_min_temp)/alpha_max)
	     end
	     #println("trial:  ", xl_trial,"  ",xu_trial)
	     xu_trial = min(xu_trial, xu)
	     xl_trial = max(xl_trial, xl)
	     
	    
	     #second round if sq > 0
	     if abs(square_coeff) > 0
	        sum_min_temp = sum_min
             	sum_max_temp = sum_max 
		sum_min_temp += min(alpha_min*xl_trial, alpha_min*xu_trial, alpha_max*xl_trial, alpha_max*xu_trial)
		sum_max_temp += max(alpha_min*xl_trial, alpha_min*xu_trial, alpha_max*xl_trial, alpha_max*xu_trial)
		sqrt_ub = sqrt(max((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff, 0))
		sqrt_lb = sqrt(max(min((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff),0))
		#println("suqare_lb", sqrt_lb, "  ", sqrt_ub)
		xu_trial = min(xu_trial, sqrt_ub)
             	xl_trial = max(xl_trial, -sqrt_ub)
		if xl_trial<= -sqrt_lb && xu_trial<= sqrt_lb && xu_trial>= -sqrt_lb
		      xu_trial = - sqrt_lb		      
		elseif xu_trial >= sqrt_lb && xl_trial<= sqrt_lb && xl_trial>= -sqrt_lb 
		      xl_trial = sqrt_lb 
		end  
		#println("square_trial:  ", xl_trial,"  ",xu_trial) 	 
	     end

	     #third round if sq > 0
             if abs(square_coeff) > 0
                sum_min_temp = sum_min
                sum_max_temp = sum_max

		#println("sum2 ", sum_min_temp, "  ",sum_max_temp)
		#println("alpha_min", alpha_min)
		#println("alpha_max", alpha_max)
		#println("square_coeff", square_coeff)
		if alpha_min <= 0 && alpha_max >=0
                   sum_min_temp += min(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4, 0)
                   sum_max_temp += max(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4, 0)
		else
		   sum_min_temp += min(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4)
                   sum_max_temp += max(-alpha_min^2/square_coeff/4, -alpha_max^2/square_coeff/4)
		end

		#println("sum", sum_min_temp, "  ",sum_max_temp)
                sqrt_ub = sqrt(max((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff, 0))
                sqrt_lb = sqrt(max(min((ub - sum_min_temp)/square_coeff, (lb - sum_max_temp)/square_coeff),0))
                #println("suqare_lb", sqrt_lb, "  ", sqrt_ub)

		temp_min = min(alpha_min/square_coeff/2, alpha_max/square_coeff/2)
		temp_max = max(alpha_min/square_coeff/2, alpha_max/square_coeff/2)
		
                xu_trial = min(xu_trial, sqrt_ub - temp_min)
                xl_trial = max(xl_trial, - sqrt_ub - temp_max)
		#x+temp_max > sqrt_lb or x+temp_min <=-sqrt_lb
		if (sqrt_lb-temp_max) > (-sqrt_lb - temp_min)
                   if xl_trial<= (-sqrt_lb - temp_min)  && xu_trial<= (sqrt_lb-temp_max) && xu_trial>= (-sqrt_lb - temp_min)
                      xu_trial = (-sqrt_lb - temp_min)
                   elseif xu_trial >= (sqrt_lb-temp_max)  && xl_trial<= (sqrt_lb-temp_max)   && xl_trial>= (-sqrt_lb - temp_min)
                      xl_trial = (sqrt_lb-temp_max)
                   end   
		end   
                #println("square_trial2:  ", xl_trial,"  ",xu_trial)
             end

	     
             if (xl_trial-xl) >= machine_error
                  setlowerbound(var, xl_trial)
                  #   println("quard   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
             end

             if (xu-xu_trial) >= machine_error
                  setupperbound(var, xu_trial)
                  #   println("quard  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
             end
          end
	  =#
	  #println("after traditional")
	  #println(P.colLower)
	  #println(P.colUpper)
	  #end
    end
    end
    
    obj = P.obj
    aff = obj.aff
    ub = U
    for j in 1:length(aff.vars)
          alpha = aff.coeffs[j]
          sum_min = aff.constant
          for k in 1:length(aff.vars)
              if k != j
                  xlk = getlowerbound(aff.vars[k])
                  xuk = getupperbound(aff.vars[k])
                  coeff = aff.coeffs[k]
                  sum_min += min(coeff*xlk, coeff*xuk)
              end
          end
          var = aff.vars[j]
          xl = getlowerbound(var)
          xu = getupperbound(var)
	  xu_trial = xu
          xl_trial = xl
          if alpha < 0
              xl_trial = (ub - sum_min)/alpha
              if (xl_trial-xl) >= machine_error
                  setlowerbound(var, xl_trial)
                  #    println("obj:   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
              end
	  elseif alpha > 0
              xu_trial = (ub - sum_min)/alpha
              if (xu-xu_trial) >= machine_error
                  setupperbound(var, xu_trial)
                  #    println("obj:  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
              end
          end
    end
    return feasible
end



function Interval_cal(aff::AffExpr)
    sum_min = aff.constant
    sum_max = aff.constant
    for k in 1:length(aff.vars)
        xlk = getlowerbound(aff.vars[k])
        xuk = getupperbound(aff.vars[k])
        coeff = aff.coeffs[k]
        sum_min += min(coeff*xlk, coeff*xuk)
        sum_max += max(coeff*xlk, coeff*xuk)
    end
    return (sum_min, sum_max)
end


function Interval_cal(quad::QuadExpr)
    qvars1 = quad.qvars1
    qvars2 = quad.qvars2
    qcoeffs = quad.qcoeffs
    aff = quad.aff

    sum_min_aff, sum_max_aff = Interval_cal(aff)
    sum_min = sum_min_aff
    sum_max = sum_max_aff
    for k in 1:length(qvars1)
        xlk1 = getlowerbound(qvars1[k])
        xuk1 = getupperbound(qvars1[k])
        xlk2 = getlowerbound(qvars2[k])
       	xuk2 = getupperbound(qvars2[k])
        coeff = qcoeffs[k]
        if qvars1[k] == qvars2[k]
            if xlk1 <= 0 && 0 <= xuk1
                sum_min += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
                sum_max += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1, 0)
            else
                sum_min += min(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
                sum_max += max(coeff*xlk1*xlk1, coeff*xuk1*xuk1)
            end
        else
            sum_min += min(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
            sum_max += max(coeff*xlk1*xlk2, coeff*xlk1*xuk2, coeff*xuk1*xlk2, coeff*xuk1*xuk2)
        end
    end
    return (sum_min, sum_max)
end


function optimality_reduction_range(P, pr, Rold, U, varsId)
    feasible = true
    #println("start medium_feasibility_reduction")
    #R = relax(P, pr, U)
    R = updaterelax(Rold, P, pr, U)
 
    left_OBBT_inner = 1		  
    left_level = 0
    #while left_level <= 0.95
	left_level = 1
        for varId in varsId
            xl = P.colLower[varId]
            xu = P.colUpper[varId]
            var = Variable(R, varId)
	    @objective(R, Min, var)
            status = solve(R)
            R_obj = getobjectivevalue(R)
            if status == :Infeasible
                feasible = false
                break
            end
            if status == :Optimal
                P.colLower[varId] = R_obj #min(R_obj, xu)
                R.colLower[varId] = R_obj #min(R_obj, xu)
            end

            R.objSense = :Max
            status = solve(R)
            R_obj = getobjectivevalue(R)
            if status == :Infeasible
                feasible = false
                break
            end	    
            if status == :Optimal
                P.colUpper[varId] = R_obj  #max(R_obj, P.colLower[varId])
                R.colUpper[varId] = R_obj  #max(R_obj, P.colLower[varId])
            end

	    if (xu - xl - P.colUpper[varId] + P.colLower[varId]) >= probing_improve
                #println("col: ", varId, "bound from [", xl," , ",xu,"] to   [", P.colLower[varId]," , ",P.colUpper[varId],"]")
		left_level = left_level * (P.colUpper[varId] - P.colLower[varId])/ (xu - xl)		
		#R=relax(P, pr, U)
                R=updaterelax(R, P, pr, U)
            end
        end	    
	left_OBBT_inner = left_OBBT_inner * left_level
	#println("optimality left_level  ", left_level, left_OBBT_inner)
    #end
    return (feasible)
end


function Sto_medium_feasibility_reduction(P, pr_children, Pex, prex, Rold, UB, LB, bVarsId)         
    feasible = true
    left_OBBT =	1
    left_OBBT_inner = 0
    #while left_OBBT_inner <= 0.1
         xlold = copy(P.colLower)
         xuold = copy(P.colUpper)

	 feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, LB, 0, true)
	 updateExtensiveBoundsFromSto!(P, Pex)

	 if feasible
	    feasible = optimality_reduction_range(Pex, prex, Rold, UB, bVarsId)	  	 
	    updateStoBoundsFromExtensive!(Pex, P) 
	 end  
	 left_OBBT_inner = 1
         for i in 1:length(P.colLower)
             if (xuold[i] + P.colLower[i] - xlold[i] - P.colUpper[i]) > small_bound_improve
                 left_OBBT_inner = left_OBBT_inner * (P.colUpper[i] - P.colLower[i])/ (xuold[i]- xlold[i])
             end
         end
	 left_OBBT = left_OBBT * left_OBBT_inner
         println("left_OBBT_all   ",left_OBBT)
    #     if !feasible
    #        break
    #     end
    #end	 	 
    return feasible
end


function Sto_slow_feasibility_reduction(P, pr_children, Pex, prex, Rold, UB, LB, bVarsId)
    feasible = true
    left_OBBT = 1
    left_OBBT_inner = 0
    while left_OBBT_inner <= 0.9
         xlold = copy(P.colLower)
         xuold = copy(P.colUpper)
         feasible = Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, UB, LB)
         println("finish fast reduction")
         updateExtensiveBoundsFromSto!(P, Pex)
         println("update fast reduction")
         if feasible
            feasible = optimality_reduction_range(Pex, prex, Rold, UB, bVarsId)
            updateStoBoundsFromExtensive!(Pex, P)
         end
         println("finish medium reduction")
         left_OBBT_inner = 1
         for i in 1:length(P.colLower)
             if (xuold[i] + P.colLower[i] - xlold[i] - P.colUpper[i]) > small_bound_improve
                 left_OBBT_inner = left_OBBT_inner * (P.colUpper[i] - P.colLower[i])/ (xuold[i]- xlold[i])
             end
         end
         left_OBBT = left_OBBT * left_OBBT_inner
         println("left_OBBT_all   ",left_OBBT)
         if !feasible
            break
         end
    end
    return feasible
end


function medium_feasibility_reduction(P, pr, U, bVarsId)
	 feasible = optimality_reduction_range(P, pr, U, bVarsId)
	  if feasible
             feasible = fast_feasibility_reduction!(P, pr, U)
         end
	 return feasible
end

function slow_feasibility_reduction(P, pr, U)
         feasible = optimality_reduction_range(P, pr, U, 1:P.numCols)
	 if feasible
	     feasible = fast_feasibility_reduction!(P, pr, U)
	 end
	 return feasible    
end


function multi_start(P, local_x, local_obj, prime_status)
      n_trial = 10      
      lb = P.colLower
      ub = P.colUpper
      ncols = length(lb)
      mlb = -1e8*ones(ncols)
      mub =  1e8*ones(ncols)
      for i in 1:ncols
          if lb[i] != -Inf
             mlb[i] = lb[i]
          end
          if ub[i] != Inf
             mub[i] = ub[i]
          end
      end
      for local_trial = 1:n_trial
          percent = rand(ncols)
          initial  = mlb + (mub-mlb).*percent
          P.colVal = initial
          prime_status_trial = solve(P)
          #println("solving Prime_trial: ",prime_status_trial)
           if prime_status_trial == :Optimal
	      local_obj_trial = getobjectivevalue(P)
              if local_obj_trial < local_obj
                  local_obj = local_obj_trial
                  prime_status = prime_status_trial
                  local_x = P.colVal
	      end	  
          end
      end
      return (local_x, local_obj, prime_status)
end


function reduced_cost_BT!(P, pr, R, U, node_L)
      mu = copy(R.redCosts)
      n_reduced_cost_BT = 0
      for varId in 1:P.numCols
          xl = P.colLower[varId]
          xu = P.colUpper[varId]
          if mu[varId] >= 1e-4
              xu_trial = xl + (U-node_L)/mu[varId]
	      #println("mu: ", mu[varId], "  col: ", varId, " [  ", xl, " , ", xu, "]", " L: ", node_L, "  U  ", U, "xu_trial ", xu_trial)
	      if (xu-xu_trial) >= machine_error
              	 P.colUpper[varId] = xu_trial
              	 R.colUpper[varId] = xu_trial
              	 #println("positive mu,  col: ", varId, "upper bound from ", xu,"  to ", xu_trial)
		 n_reduced_cost_BT += 1
              end
          elseif mu[varId] <= - 1e-4
              xl_trial = xu + (U-node_L)/mu[varId]
	      #println("mu: ", mu[varId], "  col: ", varId, " [  ", xl, " , ", xu, "]", " L: ", node_L, "  U  ", U, "xl_trial ", xl_trial)
	      if (xl_trial-xl) >= machine_error
              	 P.colLower[varId] = xl_trial
              	 R.colLower[varId] = xl_trial
              	 #println("negative mu,  col: ", varId, "lower bound from ", xl,"  to ", xl_trial)
		 n_reduced_cost_BT += 1
	      end
          end
          #=
          if (P.colLower[varId]-xl)>= large_bound_improve || (xu-P.colUpper[varId])>= large_bound_improve
              R=relax(P, pr, U)
	      #print("before add constraint in optimization based reduction")
              solve(R)
	      #println("after solve R in optimization based reduction")
              node_L = getobjectivevalue(R)
              mu = copy(R.redCosts)
          end
          =#
      end
      #feasible = fast_feasibility_reduction!(P, pr, U)
      return (n_reduced_cost_BT) #, feasible)
end

