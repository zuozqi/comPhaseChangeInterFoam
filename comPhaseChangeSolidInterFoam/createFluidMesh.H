Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);

meshPtr.reset
(
    new Foam::fvMesh
    (
        Foam::IOobject
        (
            "fluid",
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    )
);

Foam::fvMesh& mesh = meshPtr();

