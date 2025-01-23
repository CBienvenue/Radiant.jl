# Quick Start

There are many ways to install Radiant on your computer.

## Prerequisites

Before installing the Radiant package, ensure you have the following prerequisites:

1. **Julia**: Make sure Julia is installed on your system. You can download it from the [official Julia website](https://julialang.org/downloads/).

## Installation from Julia's Package Manager

This installation method is recommended for users of Radiant.

### Step 1: Open Julia REPL

First, open Julia REPL interactive session by launching Julia from a terminal or command prompt typing
```sh
> julia
```
in it.

### Step 2: Add Radiant

To add your package from the Julia prompt, you can either use
```julia
julia> using Pkg
julia> Pkg.add("Radiant")
```
or
```
julia> ]
pkg> add Radiant
```

### Step 3: Use Radiant

After installation, you can start using the Radiant package by including it in your project by adding the
```julia
using Radiant
```
line in the project's code.

## Installation from Github repository

This installation method is recommended for developers of Radiant.

### Step 1: Create a Local Repository

Create a new local repository, in which the development of new feature for Radiant will take place. Then, open the terminal or command prompt in that repository.

### Step 2: Clone the Radiant Github Repository

Clone the Radiant Github repository to your local machine by running the following command in the terminal or command prompt:
```sh
> git clone https://github.com/CBienvenue/Radiant.jl.git
```

### Step 3: Install the Radiant Package

Open Julia REPL interactive session by launching Julia from a terminal or command prompt typing
```sh
> julia
```
in it. Create a new environment in the project directory using
```
julia> ]
pkg> activate .
```
This creates a new environment in the current directory and switch to it. Your local package and its dependencies will be managed separately from the global environment. Then, navigate to the package directory by running the following command in the terminal or command prompt:
```julia
julia> cd("Radiant.jl")
```
In the Julia REPL, run
```julia
pkg> dev .
pkg> instantiate
```
This command tells Julia to develop the package located in the current directory. The dev command sets up the package for development, including adding it to your Julia environment. Then, the instantiate command is used to resolve and install any dependencies.

### Step 4: Using the Local Radiant Package

After setting up the package, you can start using it by including it in your project with the using keyword:
```julia
using Radiant
```
