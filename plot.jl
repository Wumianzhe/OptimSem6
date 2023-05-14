using Polyhedra,Plots

function prep(A, b)
  if size(b)[1] != 1
    b = permutedims(b)
  end
  A = [A[1:4, 1:2]; A[7:end, 1:2]]
  b = [b[1:4]; b[7:end]]
  return (A, b)
end

function poly(_A,_b)
    # polyhedra uses <= while my code uses >= (for now), so I need to flip signs
    local A = -permutedims(_A);
    local b = -permutedims(_b);
    local halves = [HalfSpace(A[:,i],b[i]) for i in 1:size(A)[2]]
    local p = polyhedron(hrep(halves))
    plot(p)
end
