# Ported from: http://www.math.unl.edu/~tshores1/Public/Math496S06/MatlabTools/
#                                                Examples/chap7/examp1/bvls.m
# TODO: take a critical look at speed some time
# 
# x=bvls(A,b,l,u)
#
# Solves a bounded variables least squares problem
#
#   min || Ax-b ||
#          l <= x <= u
#
# The algorithm is based on the bounded variables least squares algorithm 
# of Stark and Parker.  See
#
#  P.B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: An Algorithm
#  and Applications", Computational Statistics 10:129-141, 1995.
#

function bvls(A::AbstractMatrix, b::AbstractVector, l::AbstractVector, u::AbstractVector)
    #
    # Get the size of A for future reference.
    #
    m, n = size(A)
    #
    # Initialize a bunch of variables.  This speeds things up by avoiding 
    # realloation of storage each time an entry is added to the end of a vector.
    #
    oopslist = falses(n)
    state = fill(0, n)
    x = zeros(n)
    atbound = falses(n)
    between = falses(n)
    criti = 0
    crits = 0
    myeps = 1.0e-10
    #
    # setup an initial solution with all vars locked at a lower or upper bound.
    # set any free vars to 0.
    #
    for i in 1:n
        if ((u[i] >= +Inf) && (l[i] <= -Inf))
            x[i] = 0
            state[i] = 0
            between[i] = true
        elseif (u[i] >= +Inf)
            x[i] = l[i]
            state[i] = 1
            atbound[i] = true
        elseif (l[i] <= -Inf)
            x[i] = u[i]
            state[i] = 2
            atbound[i] = true
        elseif (abs(l[i]) <= abs(u[i]))
            x[i] = l[i]
            state[i] = 1
            atbound[i] = true
        else
            x[i] = u[i]
            state[i] = 2
            atbound[i] = true
        end
    end
    #
    # The main loop.  Stop after 10*n iterations if nothing else.  
    #
    iter = 0
    while (iter < 10*n)
        iter += 1
        #
        # optimality test.
        #
        grad = A'*(A*x-b)
        #
        # We ignore any variable that has already been tried and failed.  
        #
        grad[oopslist] = 0.
        #
        # Check for optimality.
        #
        done = true
        for i in 1:n
            if ((abs(grad[i]) > (1+norm(b))*myeps) && (state[i] == 0)) ||
                   ((grad[i] < 0) && (state[i] == 1)) ||
                   ((grad[i] > 0) && (state[i] == 2))
               done = false
               break
           end
        end 
        done && return x
        #
        # Not optimal, so we need to free up some vars at bounds.
        #
        newi = 0
        newg = 0.0
        for i in 1:length(atbound)
            if atbound[i]
                #
                # Don't free up a variable that was just locked at its bound.
                #
                i == criti && continue
                #
                # Look for the locked variable with the biggest gradient.
                #
                if (grad[i] > 0) && (state[i] == 2) && (abs(grad[i]) > newg)
                    newi = i
                    newg = abs(grad[i])
                end
                if (grad[i] < 0) && (state[i] == 1) && (abs(grad[i]) > newg)
                    newi = i
                    newg = abs(grad[i])
                end
            end
        end
        #
        # Free the locked variable with the biggest gradient if there is one.
        #
        if (newi != 0)
            atbound[newi] = false
            state[newi] = 0
            between[newi] = true
        end
        #
        # Make sure the projected problem is nontrivial.
        #
        if (length(between)==0)
          # println("Empty projected problem")
          continue
        end
        #
        # Construct the new projected problem.
        #
        Aproj = A[:, between]
        An = A[:, atbound]
        if any(atbound)
            bproj = b - An*x[atbound]
        else
            bproj = b
        end
        #
        # Solve the projected problem.
        #
        z = Aproj\bproj
        xnew = copy(x)
        xnew[between] = z

        if (newi != 0 && xnew[newi] <= l[newi] && x[newi] == l[newi]) ||
                (xnew[newi] >= u[newi] && x[newi] == u[newi])
            #
            # Ooops- the freed variable wants to go beyond its bound.  Lock it down
            # and look for some other variable.
            #
            oopslist[newi] = true
            if xnew[newi] <= l[newi] && state[newi] == 1
                state[newi] = 1
                x[newi] = l[newi]
            end
            if xnew[newi] >= u[newi] && state[newi] == 2
                state[newi] = 2
                x[newi] = u[newi]
            end
            atbound[newi] = true
            between[newi] = false
            continue
        end
        #
        # We've got a good variable freed up.  Reset the oopslist.
        #
        oopslist = falses(n)
        #
        # Move as far as possible towards the optimal solution to the projected
        # problem.
        #
        alpha = 1.
        for i in 1:length(atbound)
            if between[i]
                if (xnew[i] > u[i])
                    newalpha = min(alpha, (u[i]-x[i])/(xnew[i]-x[i]))
                    if (newalpha < alpha)
                        criti = i
                        crits = 2
                        alpha = newalpha
                    end
                end
                if (xnew[i] < l[i])
                    newalpha = min(alpha, (l[i]-x[i])/(xnew[i]-x[i]))
                    if (newalpha < alpha)
                        criti = i
                        crits = 1
                        alpha = newalpha
                    end
                end
            end
        end
        #
        # Take the step.  
        #
        x = x + alpha*(xnew-x)
        #
        # Update the state of variables.
        #
        if (alpha < 1)
            between[criti] = false
            atbound[criti] = true
            state[criti] = crits
        end

        for i in 1:n
            if (x[i] >= u[i])
                x[i] = u[i]
                state[i] = 2
                between[i] = false
                atbound[i] = true
            end
            if (x[i] <= l[i])
                x[i] = l[i]
                state[i] = 1
                between[i] = false
                atbound[i] = true
            end
        end
    end
end

