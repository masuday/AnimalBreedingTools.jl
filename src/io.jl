"""
    disp(x)
    disp(x...)
    disp(io::IO,x)
    disp(io::IO,x...)
    disp(m::Array; format)
    disp(io::IO, m::Array; format)

Write an object `x` (or objects `x...`) to `io` (or `stdout` if `io` is not given).
The output is the same as `println` and `Base.print_matrix`, but the element type will not be shown.
If an array is given, you can specify the C-style format with `format`.

```juliadoctests
A = [1 2; 3 4]
b = [9,10]

disp(A)
disp(b)
disp(A,b)
```
"""
function disp(x)
   println(x)
end
function disp(args...)
   for x in args
      disp(x)
   end
end
function disp(io::IO,x)
   println(io,x)
end
function disp(io::IO,args...)
   for x in args
      disp(io,x)
   end
end
function disp(x::Array; format::Union{Nothing,String}=nothing)
   disp(stdout,x, format=format)
end
function disp(io::IO,x::Array; format::Union{Nothing,String}=nothing)
   if isnothing(format) || ndims(x)>2
      Base.print_matrix(io,x)
      print(io,"\n")
   elseif ndims(x)==1
      for i=1:length(x)
         print_formatted(io,format,x[i])
         print("\n")
      end
   elseif ndims(x)==2
      for i=1:size(x,1)
         for j=1:size(x,2)
            print_formatted(io,format,x[i,j])
         end
         print("\n")
      end
   end
end

# see https://discourse.julialang.org/t/printf-with-variable-format-string/3805/3
print_formatted(io, fmt, args...) = @eval @printf($io, $fmt, $(args...))
