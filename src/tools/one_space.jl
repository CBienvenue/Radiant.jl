"""
    one_space(line::String) 

Transform any space in the input string into one unique space.

# Input Argument(s)
- `line::String`: string.

# Output Argument(s)
- `line_new::String`: line only with "one space" spaces.

# Reference(s)
N/A

"""
function one_space(line::String)

line_new = ""
il=1

@inbounds for i in line
    if il==1 && isspace(i)==false
        line_new = line_new * i
    else
        if isspace(i)==false || (isspace(i)==true && il>1 && isspace(line[il-1])==false)
            line_new = line_new * i
        end
    end
    il += 1
end

return line_new

end