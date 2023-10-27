#!/usr/bin/env python3

import telnetlib, socket, time, argparse

HOST = '10.0.0.241'
PORT = 4242
TIMEOUT = 9999

def send_command(cmd, host, port=PORT, timeout=TIMEOUT):
    try:
        session = telnetlib.Telnet(host, port, timeout)
    except socket.timeout:
        print("socket timeout")
    else:
        command = ('>01' + str(cmd).upper() + '\n').encode('ascii')
        session.write(command)
        time.sleep(0.5)
        response = session.read_very_eager().decode('ascii')[3:]
        #response = session.read_until(b"\n", TIMEOUT) 
        return response.strip()

def power_on(host):
    applyall = send_command('applyall', host)
    time.sleep(1)
    output = send_command('poweron', host)
    return output

def power_off(host):
    output = send_command('poweroff', host)
    return output

def get_status(host):
    status_dict = {}
    output = send_command('status', host)
    keywords = str(output).split()
    for kvpair in keywords:
        i = kvpair.split('=')
        status_dict[i[0]] = float(i[1])
    return status_dict, keywords

def get_frame(host):
    frame_dict = {}
    output = send_command('frame', host)
    keywords = str(output).split()
    for kvpair in keywords:
        i = kvpair.split('=')
        # Keys with TIME are hexadecimal according to manual
        if 'TIME' in i[0]: 
            frame_dict[i[0]] = int(i[1], 16)
        else:
            frame_dict[i[0]] = int(i[1])
    return frame_dict, keywords

def main():
    parser = argparse.ArgumentParser(description="Talk to Archon controller")
    parser.add_argument('-r', '--poweron', action='store_true', help='Power on CCD')
    parser.add_argument('-o', '--poweroff', action='store_true', help='Power off CCD')
    parser.add_argument('-s', '--status', action='store_true', help='Check status')
    parser.add_argument('-f', '--frame', action='store_true', help='Check frame')
    args = parser.parse_args()

    if args.poweron:
        print('Apply all/power on CCD')
        power_on(HOST)
    if args.poweroff:
        print('Power off CCD')
        power_off(HOST)
    if args.status:
        print('Getting status...')
        _, keywords = get_status(HOST) 
        print("\n".join(keywords))
    if args.frame:
        print('Getting frame...')
        _, keywords = get_frame(HOST) 
        print("\n".join(keywords))

if __name__ == '__main__':
    main()
