"""
    disp(x)
    disp(x...)
    disp(io::IO,x)
    disp(io::IO,x...)

Write an object `x` (or objects `x...`) to `io` (or `stdout` if `io` is not given).
The output is the same as `println` and `Base.print_matrix`, but the element type will not be shown.

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
function disp(x::Array)
   Base.print_matrix(stdout,x)
   print(stdout,"\n")
end
function disp(io::IO,x::Array)
   Base.print_matrix(io,x)
   print(io,"\n")
end
