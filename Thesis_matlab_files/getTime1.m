
function time = getTime1()

persistent h1
if isempty(h1)
    h1 = 0 ;
else
    h1 = h1 + 1 ;
end

time = h1 ;