% Spherical harmonics function integration
function intError = SphInt(quadrature, iRef, intError)

% First-order integrations
intError = FirstInt(quadrature, iRef, intError);

% Second-order integrations
intError = SecondInt(quadrature, iRef, intError);

% Third-order integrations
intError = ThirdInt(quadrature, iRef, intError);

% Fourth-order integrations
intError = FourthInt(quadrature, iRef, intError);

% Fifth-order integrations
intError = FifthInt(quadrature, iRef, intError);

% Sixth-order integrations
intError = SixthInt(quadrature, iRef, intError);

end