t:
    cargo t

alias test := t

doc:
    cargo doc --lib --no-deps

docmath:
    cargo clean --doc
    RUSTDOCFLAGS="--html-in-header doc/katex.html" cargo doc --lib --no-deps -v

changelog tag:
    git cliff -t {{tag}} -o CHANGELOG.md

release tag: (changelog tag)
    git add CHANGELOG.md
    git commit -m "chore: update changelog"
    git tag {{tag}}
    cargo publish