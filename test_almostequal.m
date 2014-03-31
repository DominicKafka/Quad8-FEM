assert(almostequal([1.01, 2, 3], [1, 2.01, 3], 0.1), ...
       'Same vectors not equal')
assert(~almostequal([2, 2, 3], [4, 2, 3], 0.1), ...
       'Dissimilar vectors compared as equal')

% test masked behaviour
assert(almostequal([1.01, 2, 3], [1, 23, 3], 0.1, [true, false, true]), ...
       'Same vectors not equal')
assert(~almostequal([2, 2, 3], [4, 2, 3], 0.1, [true, false, true]), ...
       'Dissimilar vectors compared as equal')
