
function time = getTime2()

persistent h2
if isempty(h2)
    h2 = 0 ;
else
    h2 = h2 + 1 ;
end

time = h2 ;