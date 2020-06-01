
function time = getTime3()

persistent h3
if isempty(h3)
    h3 = 0 ;
else
    h3 = h3 + 1 ;
end

time = h3 ;