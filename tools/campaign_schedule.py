"""Validate dependency barriers independently of process/resource scheduling."""


def validate_dependencies(jobs):
    names = [job['name'] for job in jobs]
    if len(names) != len(set(names)):
        raise ValueError('duplicate campaign job names')
    graph = {job['name']: job.get('depends_on', []) for job in jobs}
    for name, dependencies in graph.items():
        if not isinstance(dependencies, list) or any(d not in graph for d in dependencies):
            raise ValueError(f'invalid dependencies for {name}')
    visited, visiting = set(), set()
    def visit(name):
        if name in visiting:
            raise ValueError('cyclic campaign dependencies')
        if name in visited:
            return
        visiting.add(name)
        for dependency in graph[name]:
            visit(dependency)
        visiting.remove(name); visited.add(name)
    for name in graph:
        visit(name)


def ready(job, finished):
    complete = {j['name'] for j in finished if j.get('status') == 'complete'}
    return set(job.get('depends_on', [])).issubset(complete)
