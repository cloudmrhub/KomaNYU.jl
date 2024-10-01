# Filtering ARGS (useful for VSCode testing)
AVAILABLE_GPU_BACKENDS = ["CUDA", "AMDGPU", "Metal", "oneAPI"]
TEST_BACKENDS = filter(x->x in [AVAILABLE_GPU_BACKENDS; "CPU"], ARGS)
# Use package preferences if test_args are not specified
if isempty(TEST_BACKENDS)
    using Preferences
    TEST_BACKENDS = [load_preference(KomaNYUCore, "test_backend", "CPU")]
    @info "Using [preferences.KomaNYUCore]." backend=first(TEST_BACKENDS)
else
    @info "Using test_args" backend=first(TEST_BACKENDS)
end
# For testing with CUDA:   ] add CUDA   to KomaNYUCore/test/Project.toml
# For testing with AMDGPU: ] add AMDGPU to KomaNYUCore/test/Project.toml
# For testing with Metal:  ] add Metal  to KomaNYUCore/test/Project.toml
# For testing with oneAPI: ] add oneAPI to KomaNYUCore/test/Project.toml
USE_GPU = any(AVAILABLE_GPU_BACKENDS .∈ Ref(TEST_BACKENDS))
if "CUDA" in TEST_BACKENDS
    using CUDA
elseif "AMDGPU" in TEST_BACKENDS
    using AMDGPU
elseif "Metal" in TEST_BACKENDS
    using Metal
elseif "oneAPI" in TEST_BACKENDS
    using oneAPI
else
    @info "CPU using $(Threads.nthreads()) thread(s)" maxlog=1
end
