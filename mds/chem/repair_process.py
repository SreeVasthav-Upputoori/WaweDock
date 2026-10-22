"""Isolate native OpenMM work so the GUI can respond and cancel at any stage."""
from __future__ import annotations
import multiprocessing


def _worker(connection, plan, selected, output):
    try:
        from mds.chem.repair import build_repairs
        result=build_repairs(plan,selected,output,
            progress=lambda percent,message:connection.send(('progress',percent,message)))
        connection.send(('result',result))
    except BaseException as exc:
        connection.send(('error',f'{type(exc).__name__}: {exc}'))
    finally:
        connection.close()


def run_repair_process(plan,selected,output,ctx):
    ctx.check_cancelled()
    spawn=multiprocessing.get_context('spawn')
    receiver,sender=spawn.Pipe(duplex=False)
    process=spawn.Process(target=_worker,args=(sender,plan,list(selected),str(output)),daemon=True)
    process.start()
    sender.close()
    try:
        while True:
            ctx.check_cancelled()
            if receiver.poll(0.05):
                try: message=receiver.recv()
                except EOFError: raise RuntimeError(f'Repair worker stopped unexpectedly (exit {process.exitcode}).')
                if message[0]=='progress': ctx.progress(message[1],message[2])
                elif message[0]=='result':
                    ctx.check_cancelled()
                    return message[1]
                else: raise RuntimeError(message[1])
            elif not process.is_alive():
                raise RuntimeError(f'Repair worker exited without a result (exit {process.exitcode}).')
    finally:
        if process.is_alive(): process.terminate()
        process.join(5)
        receiver.close()
        process.close()
