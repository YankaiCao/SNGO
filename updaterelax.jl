function updaterelax(R, P, pr, U=1e10, initialValue = nothing)
    m = copyModel(R)
    #m.solver = R.solver
    m.solver = GurobiSolver(Method=2, Threads=1, Crossover=0, LogToConsole=0,OutputFlag=0, DualReductions = 0)    #m.solver=GurobiSolver(Threads=1, LogToConsole=0)
    #m.solver = GurobiSolver(Method=1, Threads=1, LogToConsole=0,OutputFlag=0)
    n = P.numCols	 
    m.colLower[1:n] = copy(P.colLower)
    m.colUpper[1:n] = copy(P.colUpper)
    #m.colVal = R.colVal[:]
    #m.colVal[1:n] = P.colVal[:]
    #m.colVal[n+1:end] = NaN
    if initialValue != nothing
       m.colVal[1:length(initialValue)] = copy(initialValue)
    end
    qbcId = m.ext[:qbcId]   #pr.qbcId 
    con = m.linconstr[length(P.linconstr)+1]    
    con.ub = U - m.obj.aff.constant

    for (key, value) in qbcId
    	    xid = key[1]
	    yid = key[2]
	    bid = value[1]
	    cid = value[2]
            xl = P.colLower[xid]
	    xu = P.colUpper[xid]
	    yl = P.colLower[yid]
	    yu = P.colUpper[yid]
	    
	    m.colLower[bid] = min(xl*yl, xl*yu, xu*yl, xu*yu)
	    m.colUpper[bid] = max(xl*yl, xl*yu, xu*yl, xu*yu)	
	    if xid == yid
	        m.colLower[bid] = max(m.colLower[bid], 0)
	    end
	    if initialValue !=nothing && length(initialValue) == n
	        m.colVal[bid] = m.colVal[xid]*m.colVal[yid]
	    end


	    if (xu !=Inf)  && (yu != Inf)
	        updateCon(m, m.linconstr[cid], - xu*yu, Inf, xid, yid, -yu, -xu)
	    end
	    if (xl !=-Inf) && (yl != -Inf)			
	        updateCon(m, m.linconstr[cid+1], - xl*yl, Inf, xid, yid, -yl, -xl)
	    end
	    if (xl !=-Inf) && (yu != Inf)			
	        updateCon(m, m.linconstr[cid+2], -Inf, -xl*yu, xid, yid, -yu, -xl)
	    end			
	    if xid != yid
	        if (xu !=Inf) && (yl != -Inf)
 	            updateCon(m, m.linconstr[cid+3], -Inf, -xu*yl, xid, yid, -yl, -xu)
		end    
	    end	    	   	    
    end
    return m
end


function updateCon(m, con, lb, ub, xid, yid, xcoeff, ycoeff)
    con.lb = lb
    con.ub = ub
    if xid == yid
       xcoeff = xcoeff + ycoeff
    end       
    definedindex = find( x->(x.col == xid), con.terms.vars)
    if length(definedindex) > 0
        con.terms.coeffs[definedindex[1]] = xcoeff
    else
        push!(con.terms.coeffs, xcoeff)
        push!(con.terms.vars, Variable(m, xid))
    end

    if xid != yid
        definedindex = find( x->(x.col == yid), con.terms.vars)
	if length(definedindex) > 0
            con.terms.coeffs[definedindex[1]] = ycoeff
        else
	    push!(con.terms.coeffs, ycoeff)
            push!(con.terms.vars, Variable(m, yid))
        end
    end
    return con
end

function addaBB!(m, pr)
    oldncon = length(m.linconstr)
    qbcId  = m.ext[:qbcId]   #pr.qbcId
    multiVariable_aBB = pr.multiVariable_aBB
    for mv in multiVariable_aBB
                terms = mv.terms
		alpha = mv.alpha
                qvars1 = terms.qvars1
                qvars2 = terms.qvars2
                qcoeffs = terms.qcoeffs
		#println("add aBB")
                newcon = LinearConstraint(AffExpr(), 0, Inf)
                aff = newcon.terms

                constant = 0
                qconstant = 0
                for k in 1:length(qvars1)
		    xid = qvars1[k].col
                    yid = qvars2[k].col
                    xv = m.colVal[xid]
                    yv = m.colVal[yid]
                    coeff = qcoeffs[k]

                    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)
                    #bid = qbcId[index[1]][3]
                    bid = mv.bilinearVars[k].col
                    bv = m.colVal[bid]

                    push!(aff.vars, Variable(m, bid))
                    push!(aff.coeffs, coeff)
                    push!(aff.vars, Variable(m,xid))
                    push!(aff.coeffs, - coeff*yv)
                    push!(aff.vars, Variable(m,yid))
                    push!(aff.coeffs, - coeff*xv)
                    constant += coeff*xv*yv
                    qconstant += coeff*bv
                end
                for k in 1:length(mv.qVars)
                    xid = mv.qVars[k].col
                    xv = m.colVal[xid]
		    bid = qbcId[xid,xid][1]
		    bv = m.colVal[bid]
		    coeff = alpha[k]

                    push!(aff.vars, Variable(m, bid))
                    push!(aff.coeffs, coeff)
                    push!(aff.vars, Variable(m,xid))
                    push!(aff.coeffs, - 2*coeff*xv)
                    constant += coeff*xv*xv
                    qconstant += coeff*bv
                end
                newcon.lb = - constant
                #if qconstant <= (constant - sigma_violation)
                   JuMP.addconstraint(m, newcon)
                #end
    end
    newncon = length(m.linconstr) - oldncon
    return newncon
end



function addOuterApproximationGrid!(m, pr, ngrid = 10)
    oldncon = length(m.linconstr)
    qbcId  = m.ext[:qbcId]   #pr.qbcId
    multiVariable_convex = pr.multiVariable_convex
    for (key, value) in qbcId
        xid = key[1]
        yid = key[2]
        bid = value[1]
        if xid == yid
            xl = m.colLower[xid]
            xu = m.colUpper[xid]
            for i = 1:ngrid
               xv = xl + (xu-xl)*i/(ngrid+1)
               @constraint(m, Variable(m, bid) - 2*xv*Variable(m, xid) + xv*xv >= 0)
            end
        end
    end
end


function addOuterApproximation!(m, pr)
    oldncon = length(m.linconstr)
    qbcId  = m.ext[:qbcId]   #pr.qbcId 
    #multiVariable_list = pr.multiVariable_list  
    multiVariable_convex = pr.multiVariable_convex

    for (key, value) in qbcId
        xid = key[1]
        yid = key[2]
        bid = value[1]
	if xid == yid
	    xv = m.colVal[xid]
	    bv = m.colVal[bid]    	    
	    if bv <= (xv^2 - sigma_violation)	    
	       @constraint(m, Variable(m, bid) - 2*xv*Variable(m, xid) + xv*xv >= 0)
	       #println("add constraint")
	    end
	end
    end

    for mv in multiVariable_convex
    #for i = 1:length(multiVariable_list)
    #    multiVariable_con = multiVariable_list[i]    
    #	for j in 1:length(multiVariable_con)
    #	    mv = multiVariable_con[j]
	    if mv.pd == 1 && length(mv.qVarsId) > 1
	        terms = mv.terms
            	qvars1 = terms.qvars1
            	qvars2 = terms.qvars2
            	qcoeffs = terms.qcoeffs
	        #println("add convex")
	        newcon = LinearConstraint(AffExpr(), 0, Inf)
                aff = newcon.terms
				
		#=					
		x_bar = m.colVal[mv.qVarsId]
		coeffs = - (mv.Q+mv.Q')*x_bar
		constant = x_bar'*coeffs/2
		constant = constant[1]
		new = [mv.bilinearVars; mv.qVars]
		aff.vars = copy([mv.bilinearVars; mv.qVars], m)
		aff.coeffs = [mv.terms.qcoeffs; coeffs]
		newcon.lb = - constant
		=#
			
		constant = 0
		qconstant = 0
	        for k in 1:length(qvars1)
		    xid = qvars1[k].col
                    yid = qvars2[k].col	  
		    xv = m.colVal[xid]
		    yv = m.colVal[yid]      
		    coeff = qcoeffs[k]
		    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)		    
		    #bid = qbcId[index[1]][3]
		    bid = mv.bilinearVars[k].col		    
		    bv = m.colVal[bid]

		    push!(aff.vars, Variable(m, bid))
		    push!(aff.coeffs, coeff)	       		    
		    push!(aff.vars, Variable(m,xid))
		    push!(aff.coeffs, - coeff*yv)
		    push!(aff.vars, Variable(m,yid))
                    push!(aff.coeffs, - coeff*xv)		  
		    constant += coeff*xv*yv
		    qconstant += coeff*bv
		end
		newcon.lb = - constant 

		if qconstant <= (constant - sigma_violation)
		   JuMP.addconstraint(m, newcon)
		end   
	    elseif mv.pd == -1 && length(mv.qVarsId) > 1
                terms = mv.terms
                qvars1 = terms.qvars1
                qvars2 = terms.qvars2
                qcoeffs = terms.qcoeffs
                #println("add convex")
                newcon = LinearConstraint(AffExpr(), -Inf, Inf)
                aff = newcon.terms

                constant = 0
                qconstant = 0
                for k in 1:length(qvars1)
                    xid = qvars1[k].col
                    yid = qvars2[k].col
                    xv = m.colVal[xid]
                    yv = m.colVal[yid]
                    coeff = qcoeffs[k]
                    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)
                    #bid = qbcId[index[1]][3]
                    bid = mv.bilinearVars[k].col
                    bv = m.colVal[bid]

                    push!(aff.vars, Variable(m, bid))
                    push!(aff.coeffs, coeff)
                    push!(aff.vars, Variable(m,xid))
                    push!(aff.coeffs, - coeff*yv)
                    push!(aff.vars, Variable(m,yid))
                    push!(aff.coeffs, - coeff*xv)
                    constant += coeff*xv*yv
                    qconstant += coeff*bv
                end
                newcon.ub = - constant

                if qconstant >= (constant + sigma_violation)
                   JuMP.addconstraint(m, newcon)
                end
	    end
#	end
    end
    newncon = length(m.linconstr) - oldncon
    return newncon
end


function checkOuterApproximation(m, pr)
    needCut = false
    qbcId  = m.ext[:qbcId] #pr.qbcId 
    for (key, value) in qbcId
        xid = key[1]
        yid = key[2]
        bid = value[1]
        cid = value[2]

        if xid == yid
            xv = m.colVal[xid]
            bv = m.colVal[bid]
            if bv <= (xv^2 - sigma_violation)
	       needCut = true
	       return needCut
            end
        end
    end

    multiVariable_convex = pr.multiVariable_convex
    #multiVariable_list = pr.multiVariable_list
    for mv in multiVariable_convex
    #for i = 1:length(multiVariable_list)
    #    multiVariable_con = multiVariable_list[i]
    #    for j in 1:length(multiVariable_con)
    #        mv = multiVariable_con[j]
            terms = mv.terms
            qvars1 = terms.qvars1
            qvars2 = terms.qvars2
            qcoeffs = terms.qcoeffs

            if mv.pd == 1 && length(mv.qVarsId) > 1
                constant = 0
                qconstant = 0
                for k in 1:length(qvars1)
                    xid = qvars1[k].col
                    yid = qvars2[k].col
                    xv = m.colVal[xid]
                    yv = m.colVal[yid]
                    coeff = qcoeffs[k]
                    bid = mv.bilinearVars[k].col
                    bv = m.colVal[bid]
                    constant += coeff*xv*yv
                    qconstant += coeff*bv
                end
                if qconstant <= (constant - sigma_violation)
		    needCut = true
               	    return needCut
                end
            elseif mv.pd == -1 && length(mv.qVarsId) > 1
                constant = 0
                qconstant = 0
                for k in 1:length(qvars1)
                    xid = qvars1[k].col
                    yid = qvars2[k].col
                    xv = m.colVal[xid]
                    yv = m.colVal[yid]
                    coeff = qcoeffs[k]
                    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)
                    #bid = qbcId[index[1]][3]
                    bid = mv.bilinearVars[k].col
                    bv = m.colVal[bid]
                    constant += coeff*xv*yv
                    qconstant += coeff*bv
                end
                if qconstant >= (constant + sigma_violation)
       		    needCut = true
                    return needCut
                end
            end
    #    end
    end
    return needCut
end

