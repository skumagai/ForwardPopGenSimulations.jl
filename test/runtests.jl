module TestForwardPopGenSimulations

using ForwardPopGenSimulations
using FactCheck

# write your own tests here
testnames = [
]

for test in tests
    Core.include("$test.jl")
end

end
