"""Accumulate nested phase counters separately from top-level scene accounting."""
from time import monotonic


class Substeps:
    def __init__(self):
        self.nodes = {}
        self.jobs = {}

    def handle(self, stats):
        event = stats.get('substep')
        if event:
            identity = event['id']
            node = self.nodes.setdefault(identity, dict(id=identity, parent=event.get('parent'),
                name=event['name'], unit=event.get('unit', 'jobs'), kind=event.get('kind', 'phase'),
                total=event.get('total'), reused=event.get('reused', 0), done=0, status='running', started=monotonic()))
            kind = event['event']
            if not kind.startswith('job_'):
                node['total'] = event.get('total')
                node['done'] = max(node['done'], int(event.get('done', 0)))
                if kind in {'end', 'failed', 'stopped'}:
                    node['status'] = 'completed' if kind == 'end' else kind
                if kind == 'advance' and event.get('job') is not None:
                    job = self.jobs.setdefault(identity + ':' + str(event['job']), {})
                    job.update(parent=identity, status='completed')
            else:
                job = self.jobs.setdefault(identity + ':' + str(event['job']), {})
                job.update(parent=identity, scene=event.get('scene', ''))
                if kind != 'job_start' or job.get('status') not in {'completed', 'finished', 'failed'}:
                    job['status'] = {'job_start': 'running', 'job_finished': 'finished', 'job_failed': 'failed'}[kind]
        if stats.get('job_id') and not event:
            job = self.jobs.setdefault(stats['job_id'], {})
            job.update(stats=dict(stats), updated=monotonic(), scene=stats.get('scene', ''),
                       parent=stats.get('parent_substep'))

    def rows(self, *, finished=False, failed=False):
        children = {}
        for node in self.nodes.values():
            parent = node['parent'] if node['parent'] in self.nodes else None
            children.setdefault(parent, []).append(node)

        def make(node, ancestors):
            if node['id'] in ancestors:
                return None
            descendants = [row for child in children.get(node['id'], [])
                           if (row := make(child, ancestors | {node['id']})) is not None]
            status = node['status']
            if status == 'running' and finished:
                status = 'failed' if failed else 'completed' if node['total'] is not None and node['done'] >= node['total'] else 'stopped'
            group = node['kind'] == 'group'
            if group and not descendants:
                return None
            run = sum(r['run'] for r in descendants) if group else int(node['total'] or 0)
            reused = sum(r['reused'] for r in descendants) if group else int(node['reused'])
            total = run + reused
            done = sum(r['done'] for r in descendants) if group else min(run, node['done'])
            own_jobs = [j for j in self.jobs.values() if j.get('parent') == node['id']]
            active = (sum(r['active'] for r in descendants) if group else
                      sum(j.get('status') == 'running' for j in own_jobs))
            if not own_jobs and not group and node['kind'] != 'jobs' and status == 'running':
                active = 1  # A serial callback bar has one coordinating operation.
            if status != 'running':
                active = 0
            if any(r['status'] == 'failed' for r in descendants) or any(j.get('status') == 'failed' for j in own_jobs):
                status = 'failed'
            elapsed = monotonic() - node['started']
            eta = (max(0, run - done) * elapsed / done if done else None) if status == 'running' else 0 if status == 'completed' else None
            counts = dict(unused=0, reused=reused, done=done, run=run, all=total)
            name = node['name']
            for prefix in sorted((self.nodes[a]['name'] for a in ancestors), key=len, reverse=True):
                if name.startswith(prefix + ' '):
                    name = name[len(prefix) + 1:]
                    break
            return dict(name=name, id=node['id'], unit=node['unit'], **counts,
                percentages={k: 100 * v / total if total else 0 for k,v in counts.items()},
                fraction_done=done / run if run else None, active=active,
                pending=node['total'] is None and not group and status == 'running',
                worker_progress=True, status=status, eta_seconds=eta, children=descendants)

        result = []
        for node in children.get(None, []):
            row = make(node, set())
            if row:
                # Function wrappers are already represented by the main step.
                result.extend(row['children'] if node['kind'] == 'group' else [row])
        return result

    def operations(self, task_id, step):
        result = []
        now = monotonic()
        for identity, job in self.jobs.items():
            if job.get('status') != 'running':
                continue
            stats = dict(job.get('stats') or dict(n=0, total=None, prefix='starting', unit='it', elapsed=0, rate=None))
            stats['elapsed'] = (stats.get('elapsed') or 0) + max(0, now - job.get('updated', now))
            rate, total = stats.get('rate'), stats.get('total')
            eta = max(0, total - stats['n']) / rate if total is not None and rate and rate > 0 else None
            result.append(dict(task_id=task_id + ':' + identity, step=step, scene=job.get('scene', ''), stats=stats, eta_seconds=eta))
        return result

    def active_workers(self):
        running = {}
        for job in self.jobs.values():
            if job.get('status') == 'running':
                parent = job.get('parent')
                running[parent] = running.get(parent, 0) + 1
        total = 0
        for identity, count in running.items():
            parent = self.nodes.get(identity, {}).get('parent')
            seen = {identity}
            while parent in self.nodes and parent not in seen:
                if running.get(parent):
                    break  # Nested work belongs to an already-counted worker.
                seen.add(parent)
                parent = self.nodes[parent]['parent']
            else:
                total += count
        return total
