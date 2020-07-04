const word solidName("solid");
Foam::autoPtr<Foam::fvMesh> solidMeshPtr(nullptr);

solidMeshPtr.reset
(
    new Foam::fvMesh
    (
        Foam::IOobject
        (
            solidName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    )
);

Foam::fvMesh& solidMesh = solidMeshPtr();

