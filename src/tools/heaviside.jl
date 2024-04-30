"""
    heaviside(x::Number)

Heaviside step function H(x).

# Input Argument(s)
- 'x::Number': input value of the function.

# Output Argument(s)
- 'H::Number': Heaviside function H(x).

# Reference(s)
N/A

"""
function heaviside(x::Number)
if (x < 0) return 0 else return 1 end
end