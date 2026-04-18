module.exports = {
    printWidth: 120,
    tabWidth: 4,
    useTabs: false,
    semi: false,
    singleQuote: false,
    trailingComma: "none",
    overrides: [
        {
            files: ["*.{json,yml,yaml}"],
            options: {
                tabWidth: 2,
                trailingComma: "none"
            }
        },
        {
            files: "*.md",
            options: {
                tabWidth: 2,
                proseWrap: "always"
            }
        }
    ]
}
