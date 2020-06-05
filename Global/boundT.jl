function Sto_fast_feasibility_reduction!(P, pr_children, Pex, prex, Rold, U, L = -1e20, ngrid = 0, solve_relax = true)
    #println("inside sto")	    
    nfirst = P.numCols
    scenarios = PlasmoOld.getchildren(P)
    nscen = length(scenarios)
    feasibility_reduced = 1e10
    rn = 1
    changed = trues(nscen)
    updateStoFirstBounds!(P)
    updateExtensiveBoundsFromSto!(P, Pex)

    while feasibility_reduced >= 1  
        feasibility_reduced = 0
        xlold = copy(P.colLower)
        xuold = copy(P.colUpper)

	if Rold != nothing
	    #R = relax(Pex, prex, U)
	    R = updaterelax(Rold, Pex, prex, U)
            if hasBin
                for i = 1:R.numCols
                    if R.colCat[i] == :Bin
                        R.colCat[i] = :Cont
                    end
                end
            end
	    
	    if ngrid > 0
	        addOuterApproximationGrid!(R, prex, ngrid)
	    end	
            fb = fast_feasibility_reduction!(R, nothing, U)
            if !fb
	        println("infeasible 1")
                return fb
            end
            Pex.colLower = R.colLower[1:length(Pex.colLower)]
            Pex.colUpper = R.colUpper[1:length(Pex.colUpper)]

	    if solve_relax && rn == 1
	        relaxed_status = solve(R)
	        relaxed_LB = getRobjective(R, relaxed_status, L)
                if relaxed_status == :Optimal && ( (U-relaxed_LB)<= mingap || (U-relaxed_LB)/abs(relaxed_LB) <= mingap)
	            #FLB = relaxed_LB
		    #println("infeasible 2")
	            return false
                elseif relaxed_status == :Infeasible
                    return false
                end
                if relaxed_status == :Optimal
                    n_reduced_cost_BT = reduced_cost_BT!(Pex, prex, R, U, relaxed_LB)
                end
	    end
	    updateStoBoundsFromExtensive!(Pex, P)
	    updateStoFirstBounds!(P)
	end

        for (idx,scenario) in enumerate(scenarios)
	    if changed[idx]
	        #println(" BT scenario     ", idx)
                fb = fast_feasibility_reduction!(scenario, pr_children[idx], 1e10)
		#println("finish fast")
            	if !fb
		    #println("infeasible 3")
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
	#println("hello 2", rn )
	
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
	#updateStoFirstBounds!(P)


	xls = Array{Float64}(nscen)
	xus = Array{Float64}(nscen)
        for j in 1:nscen
 	    xls[j] = getlowerbound(scenarios[j][:objective_value])
	    xus[j] = getupperbound(scenarios[j][:objective_value])
        end

	status = linearBackward!(xls, xus, L,  U )	
	if status == :infeasible 
	    return false
	end
	if status == :updated
            for j in 1:nscen
            var = scenarios[j][:objective_value]
            xu = getupperbound(var)
            xu_trial = xus[j]
            if (xu-xu_trial) >= small_bound_improve
                feasibility_reduced += 1
		if var.m.colCat[var.col] != :Fixed
                    setupperbound(var, xu_trial)
		end   
                changed[j] = true
                #println("obj:  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
            end
            xl = getlowerbound(var)
            xl_trial = xls[j]
            if (xl_trial-xl) >= small_bound_improve
                feasibility_reduced += 1
		if var.m.colCat[var.col] != :Fixed
                    setlowerbound(var, xl_trial)
		end   
                changed[j] = true
            end
            end
        end

        for i in 1:nfirst
	    if (xuold[i] + P.colLower[i] - xlold[i] - P.colUpper[i]) > small_bound_improve
                feasibility_reduced = feasibility_reduced + 1
            end
        end
	updateExtensiveBoundsFromSto!(P, Pex)
        rn += 1
    end
    #println("end of sto")
    return true
end





function fast_feasibility_reduction!(P, pr, U)
    n = P.numCols
    feasible = true
    if sum(P.colLower.<=P.colUpper) < n
       feasible = false
       return feasible
    end
    
    if hasBin
        for i in 1:n
            if P.colCat[i] == :Bin
	        if P.colLower[i] >= machine_error  && P.colUpper[i] == 1.0
                    fixVar(Variable(P, i), 1.0)
                elseif P.colLower[i] == 0.0  && (1.0-P.colUpper[i]) >= machine_error 
                    fixVar(Variable(P, i), 0.0)
	        elseif P.colLower[i] > machine_error  && (1.0-P.colUpper[i]) >= machine_error	
	    	    feasible = false
       		    return feasible   
                end
            end
	end
    end
    
    feasibility_reduced = 1e10
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
         #println("fast feasibiltiy_based reduction round: ", "changed bound: ", feasibility_reduced)
    end
    return feasible
end



function linearBackward!(xls, xus, lb,  ub, coeffs=nothing)
	  status = :unupdated
	  coeffs_provided = true
	  if coeffs == nothing
	     coeffs = ones(length(xls))
	     coeffs_provided = false
	  end
	  
	  hasInf = false
	  if sum(xus.>=1e6)>0 || sum(xus.<=-1e6)>0 || sum(xls.>=1e6)>0 || sum(xls.<=-1e6)>0
	      hasInf = true
	  end

	      
          if !coeffs_provided
              sum_min = sum(xls)
	      sum_max = sum(xus)
          else
	      sum_min = 0
              sum_max = 0
	      for k in 1:length(xls)
                  sum_min += min(coeffs[k]*xls[k], coeffs[k]*xus[k])
              	  sum_max += max(coeffs[k]*xls[k], coeffs[k]*xus[k])
              end
	  end
          if sum_min >= lb && sum_max <=ub
              return :unupdated
          end
	  if !hasInf
              if (sum_min-ub) >= machine_error || (lb-sum_max) >= machine_error
             	     return :infeasible
	      end
	  end   
          for j in 1:length(xls)
                  alpha = coeffs[j]
              	  xl = xls[j]
              	  xu = xus[j]
		  
		  if abs(alpha)<=1e-6
		      continue
		  end
                  xu_trial = xu
                  xl_trial = xl
		  if !hasInf
              	      if alpha > 0
                          sum_min_except = sum_min - alpha*xl
                      	  sum_max_except = sum_max - alpha*xu
                          xu_trial = (ub - sum_min_except)/alpha
                          xl_trial = (lb - sum_max_except)/alpha
              	      elseif alpha < 0
		          sum_min_except = sum_min - alpha*xu
                          sum_max_except = sum_max - alpha*xl
                          xl_trial = (ub - sum_min_except)/alpha
                          xu_trial = (lb - sum_max_except)/alpha
                      end
		  else
		      ## be careful with underflow or overflow, need to fix
		      sum_min_except = 0
		      sum_max_except = 0
		      for k in 1:length(xls)
		          if k != j
		      	      sum_min_except += min(coeffs[k]*xls[k], coeffs[k]*xus[k])		      
			      sum_max_except += max(coeffs[k]*xls[k], coeffs[k]*xus[k])
		          end
		      end
                      if alpha > 0
                          xu_trial = (ub - sum_min_except)/alpha
                          xl_trial = (lb - sum_max_except)/alpha
                      elseif alpha < 0
                          xl_trial = (ub - sum_min_except)/alpha
                          xu_trial = (lb - sum_max_except)/alpha
                      end		 
		  end 		      

                  if (xl_trial-xl) >= (machine_error)
                      xls[j] = xl_trial
                      status = :updated
                  end
                  if (xu-xu_trial) >= (machine_error)
                      xus[j] = xu_trial
                      status = :updated
                  end

		  if !hasInf
		      if (xl_trial-xl) >= (machine_error) || (xu-xu_trial) >= (machine_error)
          	          if !coeffs_provided
              	              sum_min = sum(xls)
              	 	      sum_max = sum(xus)
          	          else
              	  	      sum_min = 0
              	  	      sum_max = 0
              	  	      for k in 1:length(xls)
                     	          sum_min += min(coeffs[k]*xls[k], coeffs[k]*xus[k])
                  	          sum_max += max(coeffs[k]*xls[k], coeffs[k]*xus[k])
              		      end
		          end
          	      end
		  end
	  end
	  return status
end


function AffineBackward!(aff, lb,  ub )
    xls = getlowerbound(aff.vars)
    xus = getupperbound(aff.vars)
    coeffs = copy(aff.coeffs)
    xls = [xls; aff.constant]
    xus = [xus; aff.constant]
    coeffs = [coeffs; 1]
    status  = linearBackward!(xls, xus, lb,  ub, coeffs)

    if status == :infeasible 
       return status 
    end
    status = :unupdated      
    #if changed
        for j in 1:length(aff.vars)
            var = aff.vars[j]
	    if var.m.colCat[var.col] != :Fixed
	       if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) &&  (xls[j] - getlowerbound(var)) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xls[j] - getlowerbound(var)) >= machine_error)
                   setlowerbound(var, xls[j])
		   status = :updated
	       end
	       if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) &&  (getupperbound(var)-xus[j]) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (getupperbound(var)-xus[j]) >= machine_error)	   
                   setupperbound(var, xus[j])
		   status = :updated
	       end   
	    end
        end
    #end
    return status
end

function multiVariableForward(mv::MultiVariable)
    terms = mv.terms	
    qvars1 = terms.qvars1
    qvars2 = terms.qvars2
    qcoeffs = terms.qcoeffs
    aff = terms.aff
    
    #println("inside multivariableforward")
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


             #println("alpha_min   ",alpha_min,  "  alpha_max   ",alpha_max)
             #println("sum_min_trial  ",sum_min_trial, "  sum_max_trial   ",sum_max_trial)
	     #println("xl  ",xl, " xu   ",xu) 


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
	    #println(sum_min, "     ", sum_max)	    

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
    #println(sum_min, "     ", sum_max)


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
	    #println("a x^2 + b x   ", a,"   ",b,"    ",sum_min, "    ",sum_max)
	end
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

             if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xl_trial-xl) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xl_trial -xl) >= machine_error)
                  setlowerbound(var, xl_trial)
                  #println("quard   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
             end

	     if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xu-xu_trial) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xu-xu_trial) >= machine_error)
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

	    if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xl_trial-xl) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xl_trial -xl) >= machine_error)
                  setlowerbound(var, xl_trial)
                  #println("quard   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
            end
            if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xu-xu_trial) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xu-xu_trial) >= machine_error)
                  setupperbound(var, xu_trial)
                  #println("quard  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
            end
        end
    end
    
end


function fast_feasibility_reduction_inner!(P, pr, U)
   #println("inside fast")
    #println(P.colLower)
    #println(P.colUpper)
    feasible = true 	 
    
    if pr != nothing
        expVariable_list = pr.expVariable_list
	logVariable_list = pr.logVariable_list
	powerVariable_list = pr.powerVariable_list
	monomialVariable_list = pr.monomialVariable_list
	#=
        for (i, ev) in enumerate(expVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
	    b = ev.b
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]
            yl = P.colLower[lvarId]
            yu = P.colUpper[lvarId]

            xmin = min( exp(b*xl), exp(b*xu))
            xmax = max( exp(b*xl), exp(b*xu))

            if op == :(<=) || op == :(==)
                    if yl - xmax >= machine_error
                        return false
                    end
                    if  yu - xmax  >= small_bound_improve
                        yu = xmax + small_bound_improve
                    end
                    if  yl - xmin >= machine_error
                            #xmin = yl
			    temp = log(yl)/b          
                            if  b >= 0
                                if (temp-xl) >= small_bound_improve
                                    xl = temp - small_bound_improve
                                end
                            else
                                if (xu-temp) >= small_bound_improve
                                    xu = temp + small_bound_improve
                                end
                            end
                    end
            end
            if op == :(>=) || op == :(==)
                    if xmin - yu >=machine_error
                        return false
                    end
                    if  xmin  - yl >= small_bound_improve
                            yl = xmin - small_bound_improve
                    end
                    if  xmax - yu >= machine_error
                            #xmax = yu
			    temp = log(yu)/b
                            if b >= 0
                                if (xu-temp) >= small_bound_improve
                                    xu = temp + small_bound_improve
                                end
                            else
                                if (temp-xl) >= small_bound_improve
                                    xl = temp - small_bound_improve
                                end
                            end
                    end
            end
            P.colLower[nlvarId] = xl
            P.colUpper[nlvarId] = xu
            P.colLower[lvarId] = yl
            P.colUpper[lvarId] = yu
        end
	=#

	#lvar op log(b * nlvar)
	for (i, ev) in enumerate(logVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]
            yl = P.colLower[lvarId]
            yu = P.colUpper[lvarId]

	    #println("within BTlog    x",lvarId, "    x",nlvarId, "   ",xl, "   ",xu, "    ",yl, "   ",yu )
	    #println(b, "     ",op)

            if b >= 0 && xl <= 0
                xl = 1e-20
            elseif (b>= 0 && xu <= 0)  || (b< 0 && xl >= 0)
                #println("warning:  log(a), a cannot be negative ")
                return false
            elseif b< 0 && xu >= 0
                xu = -1e-20
            end

            xmin = min( log(b*xl), log(b*xu))
            xmax = max( log(b*xl), log(b*xu))
	    #println("xmin  xmax:   ",xmin, xmax)
	    

            if op == :(<=) || op == :(==)
                    if yl - xmax >= machine_error
                        return false
                    end
                    if  yu - xmax  >= small_bound_improve
                        yu = xmax + small_bound_improve
			#println("yu changed to ",yu)
                    end
                    if  yl - xmin >= machine_error
                            #xmin = yl
			    temp = exp(yl)/b
                            if  b >= 0
                                if (temp-xl) >= small_bound_improve
                                    xl = temp - small_bound_improve
				    #println("xl  changed to ",xl)
                                end
                            else
                                if (xu-temp) >= small_bound_improve
                                    xu = temp + small_bound_improve
				    #println("xu changed to  ", xu)
                                end
                            end
                    end
            end

            if op == :(>=) || op == :(==)
                    if xmin - yu >=machine_error
                        return false
                    end
                    if  xmin  - yl >= small_bound_improve
                            yl = xmin - small_bound_improve
			    #println("yl changed to ",yl)
                    end
                    if  xmax - yu >= machine_error
                            #xmax = yu
			    temp = exp(yu)/b
                            if b >= 0
                                if (xu-temp) >= small_bound_improve
                                    xu = temp + small_bound_improve
				    #println("xu changed to  ", xu)
                                end
                            else
                                if (temp-xl) >= small_bound_improve
                                    xl = temp - small_bound_improve
				    #println("xl changed to  ", xl)
                                end
                            end
                    end
            end
            P.colLower[nlvarId] = xl
            P.colUpper[nlvarId] = xu
            P.colLower[lvarId] = yl
            P.colUpper[lvarId] = yu
	    #println("within BT    ",nlvarId, "   ",xl, "   ",xu, "  ",yl, "  ",yu)
        end


	#=
	# lvar op (d)^(b*nlvar)
	for (i, ev) in enumerate(monomialVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
	    d = ev.d
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]
            yl = P.colLower[lvarId]
            yu = P.colUpper[lvarId]
	    
            if d <= 0 
                error("warning:  a^x, a cannot be negative, please reformulate ")
		return true
            end
	    d = d^b
            xmin = min( d^(xl), d^(xu))
            xmax = max( d^(xl), d^(xu))
	    if yl < 0
	        yl = 0
	    end
	    if yu <= 0
	        return false
	    end

            if op == :(<=) || op == :(==)
                    if yl - xmax >= machine_error
                        return false
                    end
                    if  yu - xmax  >= machine_error
                        yu = xmax
                    end
                    if  yl - xmin >= machine_error
                            #xmin = yl
                            if (d >= 1) 
                                xl = log(yl)/log(d)
                            else
                                xu = log(yl)/log(d)
                            end
                    end
            end
            if op == :(>=) || op == :(==)
                    if xmin - yu >=machine_error
                        return false
                    end
                    if  xmin  - yl >= machine_error
                            yl = xmin
                    end
                    if  xmax - yu >= machine_error
                            #xmax = yu
                            if d <= 1
                                xu = log(yu)/log(d)
                            else
                                xl = log(yu)/log(d)
                            end
                    end
            end
            P.colLower[nlvarId] = xl
            P.colUpper[nlvarId] = xu
            P.colLower[lvarId] = yl
            P.colUpper[lvarId] = yu
        end
	=#

	
	# lvar op (b*nlvar)^d
	for (i, ev) in enumerate(powerVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
	    d = ev.d
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]
            yl = P.colLower[lvarId]
            yu = P.colUpper[lvarId]

	    #println("x",lvarId, "  x",nlvarId, "   b, d:  ",b,"     ",d)
	    #d fractional # require to be positive
	    if (round(d)-d) != 0  
		if b >= 0 && xl <= 0
		   xl = 0
		elseif (b>= 0 && xu <= 0)  || (b< 0 && xl >= 0)
		   #println("warning:  (a)^(fractional), a cannot be negative ")
		   return false
		elseif b< 0 && xu >= 0
		   xu = 0
		end   
	    end			
	    #println("b  xl  xu", b, "    ",xl, "   ",xu, "   ", yl, "   ",yu)	    
	

	    if d == 0
	        error("a^0 is detected, reformulate!")
	    end

	    if  xl<=0 && xu >= 0 && d<0
		 #println("warning: a^(-) for a = 0 is possible")
            end


	    #if ((round(d)-d) != 0) || (d%2 == 1)  || (d%2 == 0 && b*xl >=0)		
	    if !(xl <0 && xu > 0 && d<0 && d%2 != 0)	      
		if bounded(xl) && bounded(xu)
                    xmin = min( (b*xl)^d, (b*xu)^d)
                    xmax = max( (b*xl)^d, (b*xu)^d)		
		elseif bounded(xl)
		    xmin = min( (b*xl)^d, (b*Inf)^d)
		    xmax = max( (b*xl)^d, (b*Inf)^d)    		       
                elseif bounded(xu)
		    xmin = min( (-b*Inf)^d, (b*xu)^d)
		    xmax = max( (-b*Inf)^d, (b*xu)^d)
		else
		    xmin = -Inf
		    xmax = Inf
		end    
		
		if (d%2 == 0) && xl<=0 && xu >= 0 && d>0
		    xmin = min(xmin, 0)
		end
		if (xl <0 && xu > 0 && d<0)
		    xmax = Inf
		end
		#println("xmin  xmax    ", xmin, "    ",xmax)
				

		if op == :(<=) || op == :(==)
                    if yl - xmax >= machine_error 
                        return false
                    end
                    if  yu - xmax  >= small_bound_improve #machine_error
                            yu = xmax + small_bound_improve
			    #println("xmax   ", xmax," yu changed to  ", yu )
                    end
                    if  yl - xmin >= machine_error && bounded(yl)
		    	    #xmin = yl			    
			    if (d%2 != 0)
			        temp = yl^(1/d) /b
                                if  ((b >= 0 && d >=0) || (b<=0 && d<=0)) 
				    if (temp-xl) >= small_bound_improve
                                        xl = temp - small_bound_improve
				    end	
                                else
				    if (xu-temp) >= small_bound_improve
                                        xu = temp + small_bound_improve
				    end	
				end
			    else
				temp  = yl^(1/d)/abs(b)
			        if d >= 0
				     # x cannot be within -temp; temp				    			    
				     if xu >= temp && xl >= -temp && (temp-xl) >= small_bound_improve
				         xl = temp - small_bound_improve
				     end
				     if xl <= -temp  && xu<= temp  && (xu+temp) >= small_bound_improve
				         xu = -temp + small_bound_improve
				     end	 
				else
				     # x must be within -temp ; temp	
				     if  -temp-xl >= small_bound_improve
				     	 xl = -temp - small_bound_improve
				     end	 	
				     if xu - temp >= small_bound_improve
				     	 xu = temp + small_bound_improve				     	 
				     end	 
                                end
			    end	
                    end
		end    

                if op == :(>=) || op == :(==)		   
                    if xmin - yu >= machine_error
		        #println("xmin   ", xmin," yu  ", yu, "infeasible " )
                        return false
                    end
                    if  xmin  - yl >= small_bound_improve #machine_error
                            yl = xmin - small_bound_improve
			    #println("xmin   ", xmin," yl  ", yl )
                    end
		    
                    if  xmax - yu >= machine_error  && bounded(yu)
		    	    #xmax = yu
			    if (d%2 != 0)
                                if  ( ((b >= 0 && d >=0) || (b<=0 && d<=0)))				
                                    xu_trial = yu^(1/d)/b
				    if (xu-xu_trial) >= small_bound_improve
				        xu = xu_trial + small_bound_improve
				    end				
                                else
				    xl_trial = yu^(1/d)/b
				    if (xl_trial-xl) >= small_bound_improve	
                                        xl = xl_trial - small_bound_improve
				    end	
			        end
			    else
				temp = yu^(1/d)/abs(b)
		    		if d >= 0
				    # x must be within -temp ; temp
                                    if -xl - temp >=  small_bound_improve
                                        xl  = - temp - small_bound_improve
                                    end
                                    if  xu - temp >= small_bound_improve
                                        xu = temp + small_bound_improve
                                    end
                                else
                                     # x cannot be within -temp; temp
                                     if xu >= temp && xl >= -temp &&  (temp-xl) >= small_bound_improve 
                                         xl = temp - small_bound_improve
                                     end
                                     if xl <= -temp && xu<= temp  &&  (xu+temp) >= small_bound_improve
                                         xu = -temp + small_bound_improve
                                     end
                                end	
                            end
                    end		    
		end    	
	    elseif  xl<=0 && xu >= 0 && d<0 && d%2 == 0
                xmin = min( (b*xl)^d, (b*xu)^d)
                if op == :(<=) || op == :(==)
                    if  yl - xmin >= machine_error
                            #xmin = yl
			    temp  = yl^(1/d)/abs(b)	  
			    if xu - temp >= small_bound_improve
                                xu = temp + small_bound_improve
			    end
			    if - xl - temp >= - small_bound_improve
			        xl = -temp - small_bound_improve
                            end
                    end
                elseif op == :(>=) || op == :(==)
                    if xmin - yu >=machine_error
                        return false
                    end
                    if  xmin  - yl >= machine_error
                            yl = xmin
                    end
                end
            end
	    #println("b  xl  xu", b, "    ",xl, "   ",xu, "   ", yl, "   ",yu)	  
            P.colLower[nlvarId] = xl
            P.colUpper[nlvarId] = xu
            P.colLower[lvarId] = yl
            P.colUpper[lvarId] = yu
	end	

    end
    

    for i = 1:length(P.linconstr)
      	  con = P.linconstr[i] 
	  aff = con.terms

	  status = AffineBackward!(aff,  con.lb,  con.ub)	 
	  if status == :infeasible 
             feasible = false
             #println("infeasible   ", con)
             return (feasible)	  
	  end   
    end
    

    if length(P.quadconstr) > 0
    multiVariable_list = pr.multiVariable_list 

    #println("working on quad")
    #println(length(P.quadconstr))
    #println(multiVariable_list)
   

    for i = 1:length(P.quadconstr)
          con = P.quadconstr[i]
          terms = con.terms
          qvars1 = terms.qvars1
          qvars2 = terms.qvars2
          qcoeffs = terms.qcoeffs
          aff = terms.aff

	  #=
	  println(con, i)
	  println("mvs1    ",multiVariable_list[i])
	  for k in 1:length(qcoeffs)
	      println("Qid    ", qvars1[k].col, "   ",qvars2[k].col, "   ", qcoeffs[k])
	      println("x1  lb ub  ",P.colLower[qvars1[k].col], "    ", P.colUpper[qvars1[k].col]) 
	      println("x2  lb ub  ",P.colLower[qvars2[k].col], "    ", P.colUpper[qvars2[k].col])
	  end
	  for k in 1:length(aff.vars)
	      println("aid   ", aff.vars[k].col, "   ", aff.coeffs[k])
	      println("x  lb   ub ", P.colLower[aff.vars[k].col], "    ", P.colUpper[aff.vars[k].col])
	  end
	  =#

          lb = 0
          ub = 0
	  if con.sense == :(<=)
	     lb = -1e20
	  end
	  if con.sense == :(>=)
	     ub = 1e20
	  end
	  #println(terms, "   ", lb, "   ",ub)	  	  


	  mv_con = copy(multiVariable_list[i], P)
	  mvs =	mv_con.mvs
	  remainaff = mv_con.aff
	  nmw = length(mv_con.mvs)

	  #=
	  println(mvs)
	  println("remainaff",remainaff)
	  =#

	  mv_sum_min = Array{Float64}(nmw + 1)
	  mv_sum_max = Array{Float64}(nmw + 1)	  
	  for j = 1:nmw
	      mv_sum_min[j], mv_sum_max[j] = multiVariableForward(mvs[j])
	      #println("result of multi forward  ", j, "  ",mv_sum_min[j], "    ",mv_sum_max[j] )
	  end
	  mv_sum_min[nmw+1], mv_sum_max[nmw+1] = Interval_cal(remainaff)
	  #println("result of interval cal:  ", mv_sum_min[nmw+1], "    " ,mv_sum_max[nmw+1])


	  mv_lb = copy(mv_sum_min)
	  mv_ub = copy(mv_sum_max)

          status = linearBackward!(mv_lb, mv_ub, lb,  ub)

	  # check if constraint is feasible
          if status == :infeasible
             feasible = false
             return (feasible)
          end

	  #if status == :updatedd	  
          for j = 1:nmw
              multiVariableBackward!(mvs[j],mv_sum_min[j], mv_sum_max[j], mv_lb[j], mv_ub[j])
          end

	  AffineBackward!(remainaff, mv_lb[end],  mv_ub[end])


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

            if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xl_trial-xl) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xl_trial -xl)>= machine_error)
                  setlowerbound(var, xl_trial)
                  #   println("quard   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
             end


            if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xu-xu_trial) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xu-xu_trial) >= machine_error)
                  setupperbound(var, xu_trial)
                  #   println("quard  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
             end
          end
	  =#
    end
    end
    
    #println("obj")
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
  	      if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xl_trial-xl) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xl_trial -xl)>= machine_error)
                  setlowerbound(var, xl_trial)
                  #    println("obj:   col: ", var.col,  " lower bound from ", xl,"   to   ",xl_trial)
              end
	  elseif alpha > 0
              xu_trial = (ub - sum_min)/alpha

	      if ((var.m.colCat[var.col] == :Bin || var.m.colCat[var.col] == :Int ) && (xu-xu_trial) >= small_bound_improve) || (var.m.colCat[var.col] == :Cont && (xu-xu_trial) >= machine_error)
                  setupperbound(var, xu_trial)
                  #    println("obj:  col: ", var.col,  " upper bound from ", xu,"   to   ",xu_trial)
              end
          end
    end
    #println("end of obj")
    return feasible
end



function Interval_cal(aff::AffExpr)
    sum_min = aff.constant
    sum_max = aff.constant
    #println("inside interval_cal:      ", aff)

    for k in 1:length(aff.vars)
        xlk = getlowerbound(aff.vars[k])
        xuk = getupperbound(aff.vars[k])
        coeff = aff.coeffs[k]
        sum_min += min(coeff*xlk, coeff*xuk)
        sum_max += max(coeff*xlk, coeff*xuk)
	#println("interval_cal   ", aff.vars[k].col, "   xl    ", xlk, "  xu  ", xuk, "  coeff  ", coeff, "   sum_min  ",sum_min, "    ", sum_max)
    end
    return (sum_min, sum_max)
end


function Interval_cal(quad::QuadExpr)
    #println("inside interval cal")
    qvars1 = quad.qvars1
    qvars2 = quad.qvars2
    qcoeffs = quad.qcoeffs
    aff = quad.aff

    sum_min_aff, sum_max_aff = Interval_cal(aff)
    sum_min = sum_min_aff
    sum_max = sum_max_aff
    #println("sum min  max:   ", sum_min, "   " ,sum_max )
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
	#println("xl xu xl2 xu2   ", xlk1, "  ",xuk1,"  ",xlk2, "  ",xuk2)
	#println("sum min  max:   ", sum_min, "   " ,sum_max )
    end
    return (sum_min, sum_max)
end


function optimality_reduction_range(P, pr, Rold, U, varsId)
    if !OBBT
        return true
    end
    feasible = true
    #println("start medium_feasibility_reduction")
    #R = relax(P, pr, U)
    R = updaterelax(Rold, P, pr, U)

    if hasBin
                for i = 1:R.numCols
                    if R.colCat[i] == :Bin
                        R.colCat[i] = :Cont
                    end
                end
    end
 
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
	        if ( (P.colCat[varId] == :Bin || P.colCat[varId]== :Int ) &&   (R_obj-xl)>= small_bound_improve) || (P.colCat[varId] == :Cont &&  (R_obj-xl)>= machine_error)
                    P.colLower[varId] = R_obj #min(R_obj, xu)
                    R.colLower[varId] = R_obj #min(R_obj, xu)
		end
            end

            R.objSense = :Max
            status = solve(R)
            R_obj = getobjectivevalue(R)
            if status == :Infeasible
                feasible = false
                break
            end	    

            if status == :Optimal 
	        if  (  (P.colCat[varId] == :Bin || P.colCat[varId]== :Int ) &&   (xu-R_obj)>= small_bound_improve) || (P.colCat[varId] == :Cont &&  (xu-R_obj)>= machine_error)
                    P.colUpper[varId] = R_obj  #max(R_obj, P.colLower[varId])
                    R.colUpper[varId] = R_obj  #max(R_obj, P.colLower[varId])		
		end    
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
         #println("left_OBBT_all   ",left_OBBT)
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
         #println("finish fast reduction")
         updateExtensiveBoundsFromSto!(P, Pex)
         #println("update fast reduction")
         if feasible
            feasible = optimality_reduction_range(Pex, prex, Rold, UB, bVarsId)
            updateStoBoundsFromExtensive!(Pex, P)
         end
         #println("finish medium reduction")
         left_OBBT_inner = 1
         for i in 1:length(P.colLower)
             if (xuold[i] + P.colLower[i] - xlold[i] - P.colUpper[i]) > small_bound_improve
                 left_OBBT_inner = left_OBBT_inner * (P.colUpper[i] - P.colLower[i])/ (xuold[i]- xlold[i])
             end
         end
         left_OBBT = left_OBBT * left_OBBT_inner
         #println("left_OBBT_all   ",left_OBBT)
         if !feasible
            break
         end
    end
    return feasible
end



function reduced_cost_BT!(P, pr, R, U, node_L)
      mu = copy(R.redCosts)
      n_reduced_cost_BT = 0
      if length(mu) == 0 
      	  return n_reduced_cost_BT
      end
      
      for varId in 1:P.numCols
          xl = P.colLower[varId]
          xu = P.colUpper[varId]
          if mu[varId] >= 1e-4
              xu_trial = xl + (U-node_L)/mu[varId]
	      #println("mu: ", mu[varId], "  col: ", varId, " [  ", xl, " , ", xu, "]", " L: ", node_L, "  U  ", U, "xu_trial ", xu_trial)

	      
	      if ((P.colCat[varId] == :Bin || P.colCat[varId]== :Int) &&  (xu-xu_trial) >= small_bound_improve )|| (P.colCat[varId] == :Cont && (xu-xu_trial) >= machine_error)
              	 P.colUpper[varId] = xu_trial
              	 R.colUpper[varId] = xu_trial
              	 #println("positive mu,  col: ", varId, "upper bound from ", xu,"  to ", xu_trial)
		 n_reduced_cost_BT += 1
              end
          elseif mu[varId] <= - 1e-4
              xl_trial = xu + (U-node_L)/mu[varId]
	      #println("mu: ", mu[varId], "  col: ", varId, " [  ", xl, " , ", xu, "]", " L: ", node_L, "  U  ", U, "xl_trial ", xl_trial)
	      if ((P.colCat[varId] == :Bin || P.colCat[varId]== :Int) && (xl_trial-xl) >= small_bound_improve) || (P.colCat[varId] == :Cont && (xl_trial-xl) >= machine_error)
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

