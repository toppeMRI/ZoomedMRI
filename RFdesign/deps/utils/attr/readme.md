# An Unified Interface of Parsing and Set/Getting for Matlab

`attrParser.m`,\
an alternative to `vararg_pair.m`. Core part based on MatLab's built-in parsing,
which might have better performance as Matlab evolves overtime.\
`attrs.m`,\
an unified interface for `fieldnames()` and `properties()`. It retrieves fields
regardless of whether its input is a structure or an object.\
`chkattrs.m`,\
assert whether a field exists in its input.\
`cutattrs.m`,\
cut prescribed fields from its input.\
`getattrs.m`,\
retrieve (multiple) fields values from its input, e.g.,
`[fa, fb] = getattrs(input, {'fa','fb'});`\
`isattr.m`,\
an unified interface for `isfield()` and `isprop()`.\
`mrgattrs.m`,\
merge fields in its src_intpu to its dst_input.\
`rmattrs.m`,\
remove prescribed fields from its input.