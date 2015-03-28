module TestForwardPopGenSimulations

using ForwardPopGenSimulations
using FactCheck

# write your own tests here
testnames = [
    "chromosomalstorages"
]

for test in testnames
    Core.include("$test.jl")
end

end
