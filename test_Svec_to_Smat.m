function test_Svec_to_Smat()
Svec = 1:3;

Smat = Svec_to_Smat(Svec);

% Smat must be symmetrical
assert(all(all(Smat == Smat')), ...
       'Smat is not symmetrical');
% Smat must be 4x4
assert(all(size(Smat) == [4, 4]), ...
       'Smat must be 4x4')