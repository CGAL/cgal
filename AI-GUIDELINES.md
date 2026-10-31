# Guidelines for AI Usage in the CGAL Project

This document provides guidance for AI tools and developers using AI
assistance when contributing to the CGAL project.

## General Guidelines

 - ALL code, including AI-generated or assisted code, MUST comply with our code formatting conventions, code naming conventions, and commit guidelines.

 - Know your personal limits.
   AI assistance is a tool helping with development, not a replacement for development skills.

 - Only use AI assistance in parts of the codebase you are familiar with.
   You MUST test your changes to the best of your ability BEFORE submitting your code.

 - No slop. Do not submit code changes that solely rely on "vibe coding" without
   understanding the generated code.

 - AI must be restrained to write the appropriate amount of code comments. Lengthy LLM-oriented explanations are forbidden.

 - AI assistance MUST be disclosed in the commit message. See "Attribution".
   Using AI assistance repeatedly without disclosure or with the goal of appearing
   more capable than you actually are WILL have consequences, up to being reported.

 - Use of AI-assisted tooling to write extensive parts of code is reserved for established team members who are experienced with the codebase and its conventions.

 - Every commit and PR must be done by humans and with proper scoping. All messages must be written by a human being.

 - Code MUST NOT be authored by an AI agent. Even with using AI assistance, you are the sole author of the code you submit. AI agents must never have (Co-) authorship of your code.

## Licensing and Legal Requirements

All contributions must comply with our licensing requirements.
This means that all code must be compatible with the original license of the package/file (GPLv3+, LGPLv3+, ...)

## Developer's Responsibilities

AI agents cannot be used as an excuse for low quality. The human submitter
is responsible for:

- Reviewing all AI-generated code
- Ensuring compliance with licensing requirements
- Taking full responsibility for the contribution

## Attribution

When AI tools contribute to the project’s development, proper attribution
helps track the evolving role of AI in the development process.
Contributions should include an Assisted-by tag in the following format:

```
Assisted-by: AGENT_NAME:MODEL_VERSION [TOOL1] [TOOL2]
```

Where:

- `AGENT_NAME` is the name of the AI tool or framework
- `MODEL_VERSION` is the specific model version used
- `[TOOL1] [TOOL2]` are optional specialized analysis tools used
  (e.g., coccinelle, sparse, smatch, clang-tidy)

Basic development tools (git, gcc, make, editors) should not be listed.

Example:

```
Assisted-by: Claude:claude-3-opus coccinelle sparse [OR]
Assisted-by: Claude:claude-opus-4.6

```

