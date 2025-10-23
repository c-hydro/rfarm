# CODEOWNERS

This file defines the owners responsible for specific parts of the repository.

For more information, see:  
🔗 [GitHub Docs – About CODEOWNERS](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-code-owners#example-of-a-codeowners-file)

---

## Repository Ownership

You can assign individuals or teams as code owners.  
Code owners automatically receive review requests when someone opens a pull request that modifies files they own.

### Example

```text
# Each line is a file pattern followed by one or more GitHub usernames or team names.
# The following assigns ownership of all files in the repository:
* @fabiodelogu
```

---

### Notes

- Code owners must have **write** permissions in the repository.  
- Patterns follow the same syntax used in `.gitignore` files.  
- You can define ownership per folder or file as needed.
