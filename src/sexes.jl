export SexType,

       Female,
       Male,
       Mother,
       Father,
       Daughter,
       Son,

       value

immutable SexType
    data::Int
end

@inline value(s::SexType) = s.data

const Female = SexType(1)
const Male = SexType(2)
# Mother, Father, Daughter, and Son are alias to Female and Male. These aliases are
# used in parameter specification to reduce confusion.
const Mother = SexType(1)
const Father = SexType(2)
const Daughter = SexType(1)
const Son = SexType(2)
