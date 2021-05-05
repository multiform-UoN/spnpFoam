namespace electrokineticConstants
{
    const dimensionedScalar epsilonZero_
    (
            "vacuumPermittivity",
            dimensionSet(-1, -3, 4, 0, 0, 2, 0),
            8.8541878176e-12
    );


    const dimensionedScalar Na_
    (
            "AvogradoNumber",
            dimensionSet(0, 0, 0, 0, -1, 0, 0),
            6.022140857e+23
    );

    const dimensionedScalar e_
    (
            "elementaryCharge",
            dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),
            1.6021766208e-19
    );

    const dimensionedScalar k_
    (
            "BoltzmannConstant",
            dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),
            1.38064852e-23
    );

    const dimensionedScalar Fy_
    (
           "FaradayConstant",
           (e_ * Na_).dimensions(),
           (e_ * Na_).value()
    );
}
